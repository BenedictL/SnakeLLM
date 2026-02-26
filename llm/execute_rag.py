"""
llm/execute_rag.py
==================
EXECUTE RAG — Stage 2 retrieval
"What is the exact container URI, CLI syntax, and resource requirement for this tool?"

Retrieves precise, version-pinned execution details from the BioContainers index.
This grounds the LLM in real container URIs — preventing hallucinated image names.

Knowledge base: data/biocontainers/*.json  (built by biocontainers_indexer.py)
"""

from __future__ import annotations
import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
from rank_bm25 import BM25Okapi
import chromadb
from sentence_transformers import SentenceTransformer

log = logging.getLogger(__name__)

BIOCONTAINERS_DIR = Path("data/biocontainers")
CHROMA_DIR        = Path("data/chroma/execute_rag")
EMBED_MODEL       = "all-MiniLM-L6-v2"   # fast; good for short tool descriptions


class ExecuteRAG:
    """
    Retrieves exact execution details for specific tools:
    - Pinned container URIs  (quay.io/biocontainers/star:2.7.10a--h9ee0642_1)
    - Verified CLI flags and syntax
    - Resource benchmarks (CPU, RAM, typical runtime)

    Two retrieval paths:
    1. Exact match lookup  — if the user names a tool directly (e.g. "STAR")
    2. Semantic fallback   — if the tool is described but not named
    """

    def __init__(self, embed_model: str = EMBED_MODEL):
        log.info(f"Initializing Execute RAG (embed model: {embed_model})")
        self.embed_model = SentenceTransformer(embed_model)

        # Load all indexed BioContainers records into memory for fast lookup
        self.tool_registry: dict[str, dict] = self._load_tool_registry()

        self._load_chromadb()
        self._load_bm25_index()
        log.info(f"  ✓ Execute RAG ready ({len(self.tool_registry)} tools in registry)")

    # ── INITIALIZATION ────────────────────────────────────────────────────────

    def _load_tool_registry(self) -> dict[str, dict]:
        """
        Loads all BioContainers JSON records from disk into a dict keyed by tool name.
        Also builds name aliases (e.g. "bioconductor-deseq2" → "deseq2").
        """
        registry = {}
        if not BIOCONTAINERS_DIR.exists():
            log.warning(f"BioContainers index not found at {BIOCONTAINERS_DIR}. "
                        f"Run: python -m llm.biocontainers_indexer --defaults")
            return registry

        for f in BIOCONTAINERS_DIR.glob("*.json"):
            try:
                record = json.loads(f.read_text())
                name   = record["tool_name"].lower()
                registry[name] = record

                # Build aliases for easier lookup
                # "bioconductor-deseq2" → also accessible as "deseq2"
                if "-" in name:
                    parts = name.split("-")
                    registry[parts[-1]] = record   # last segment as alias
            except Exception as e:
                log.warning(f"  Could not load {f}: {e}")

        # Hard-coded aliases for tools stored under a different name
        if "gatk4" in registry:
            registry["gatk"] = registry["gatk4"]
            
        if "subread" in registry:
            registry["featurecounts"] = registry["subread"]
            registry["featureCounts"] = registry["subread"]

        return registry

    def _load_chromadb(self):
        """Load ChromaDB vector store for semantic fallback search."""
        CHROMA_DIR.mkdir(parents=True, exist_ok=True)
        self.chroma_client = chromadb.PersistentClient(path=str(CHROMA_DIR))
        self.collection    = self.chroma_client.get_or_create_collection(
            name="execute_rag",
            metadata={"hnsw:space": "cosine"}
        )
        count = self.collection.count()
        log.info(f"  ✓ Execute RAG ChromaDB loaded ({count} documents indexed)")

    def _load_bm25_index(self):
        """Build BM25 index over tool embedding texts for keyword search."""
        self.bm25_docs     = []
        self.bm25_tool_ids = []

        for name, record in self.tool_registry.items():
            text = record.get("embedding_text", record.get("description", name))
            self.bm25_docs.append(text)
            self.bm25_tool_ids.append(name)

        if self.bm25_docs:
            tokenized = [doc.lower().split() for doc in self.bm25_docs]
            self.bm25 = BM25Okapi(tokenized)
        else:
            self.bm25 = None

    # ── INDEXING ──────────────────────────────────────────────────────────────

    def index_documents(self, force_reindex: bool = False):
        """Index all BioContainers records into ChromaDB for semantic search."""
        if self.collection.count() > 0 and not force_reindex:
            log.info("Execute RAG already indexed. Pass force_reindex=True to rebuild.")
            return

        if not self.tool_registry:
            log.error("No tools in registry. Run biocontainers_indexer.py first.")
            return

        log.info(f"Indexing {len(self.tool_registry)} tools into Execute RAG...")

        docs, ids, metas, embeddings = [], [], [], []
        seen_ids = set()

        for name, record in self.tool_registry.items():
            # Skip duplicates caused by aliases
            if name in seen_ids:
                continue
            seen_ids.add(name)

            text = record.get("embedding_text", record.get("description", name))
            emb  = self.embed_model.encode(text).tolist()

            docs.append(text)
            ids.append(f"tool::{name}")
            metas.append({
                "tool_name":    record["tool_name"],
                "display_name": record.get("display_name", name),
                "container_uri": record.get("best_container", {}).get("full_uri", ""),
            })
            embeddings.append(emb)

        self.collection.upsert(
            documents=docs,
            embeddings=embeddings,
            ids=ids,
            metadatas=metas,
        )
        log.info(f"  ✓ Indexed {len(docs)} tools into Execute RAG")

    # ── RETRIEVAL ─────────────────────────────────────────────────────────────

    def retrieve_for_tools(self, tool_names: list[str]) -> list[dict]:
        """
        Primary path: retrieve execution details for a known list of tool names.
        Uses exact match first, then semantic fallback for each tool.
        """
        results = []
        for name in tool_names:
            record = self._exact_lookup(name)
            if record:
                results.append({"tool": name, "record": record, "method": "exact"})
            else:
                # Semantic fallback
                hits = self._semantic_search(name, top_k=1)
                if hits:
                    results.append({"tool": name, "record": hits[0], "method": "semantic_fallback"})
                else:
                    log.warning(f"  No container found for tool: '{name}'")
                    results.append({"tool": name, "record": None, "method": "not_found"})

        return results

    def retrieve_by_description(self, description: str, top_k: int = 5) -> list[dict]:
        """
        Fallback path: find tools by semantic description when tool name is unknown.
        e.g. "fast short read aligner for DNA" → BWA-MEM2
        """
        semantic = self._semantic_search(description, top_k=top_k)
        bm25     = self._bm25_search(description, top_k=top_k)
        fused    = self._rrf_fusion(semantic, bm25, final_top_k=top_k)
        return fused

    def _exact_lookup(self, tool_name: str) -> Optional[dict]:
        """Fast O(1) lookup by tool name (case-insensitive, strips whitespace)."""
        key = tool_name.lower().strip()
        return self.tool_registry.get(key, None)

    def _semantic_search(self, query: str, top_k: int) -> list[dict]:
        """Dense vector search for tool descriptions."""
        if self.collection.count() == 0:
            return []

        emb     = self.embed_model.encode(query).tolist()
        results = self.collection.query(
            query_embeddings=[emb],
            n_results=min(top_k, self.collection.count()),
            include=["documents", "metadatas", "distances"]
        )

        hits = []
        for i in range(len(results["documents"][0])):
            tool_name = results["metadatas"][0][i]["tool_name"]
            record    = self.tool_registry.get(tool_name)
            if record:
                hits.append({**record, "retrieval_score": 1 - results["distances"][0][i]})

        return hits

    def _bm25_search(self, query: str, top_k: int) -> list[dict]:
        """BM25 keyword search over tool descriptions."""
        if self.bm25 is None:
            return []

        tokens  = query.lower().split()
        scores  = self.bm25.get_scores(tokens)
        top_idx = np.argsort(scores)[::-1][:top_k]

        return [
            {**self.tool_registry[self.bm25_tool_ids[i]], "retrieval_score": float(scores[i])}
            for i in top_idx
            if scores[i] > 0 and self.bm25_tool_ids[i] in self.tool_registry
        ]

    def _rrf_fusion(
        self,
        semantic_hits: list[dict],
        bm25_hits: list[dict],
        k: int = 60,
        final_top_k: int = 5
    ) -> list[dict]:
        """Reciprocal Rank Fusion over two result lists."""
        scores    = {}
        doc_store = {}

        for rank, hit in enumerate(semantic_hits):
            tool_id = hit.get("tool_name", f"sem_{rank}")
            scores[tool_id]    = scores.get(tool_id, 0) + 1 / (k + rank + 1)
            doc_store[tool_id] = hit

        for rank, hit in enumerate(bm25_hits):
            tool_id = hit.get("tool_name", f"bm25_{rank}")
            scores[tool_id]    = scores.get(tool_id, 0) + 1 / (k + rank + 1)
            doc_store[tool_id] = hit

        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        return [doc_store[tid] for tid, _ in ranked[:final_top_k]]

    # ── CONTAINER LOOKUP ─────────────────────────────────────────────────────

    def get_container_uri(self, tool_name: str) -> str:
        """
        Direct container URI lookup for a tool name.
        Returns the best available container URI, or a placeholder if not found.
        """
        record = self._exact_lookup(tool_name)
        if record:
            uri = record.get("best_container", {}).get("full_uri", "")
            if uri:
                return uri

        log.warning(f"No container URI found for '{tool_name}'. Using placeholder.")
        return f"quay.io/biocontainers/{tool_name.lower()}:latest"

    # ── CONTEXT FORMATTING ────────────────────────────────────────────────────

    def format_for_prompt(self, tool_results: list[dict]) -> str:
        """
        Formats Execute RAG results for injection into LLM prompt.
        Focuses on container URIs and verified tool details.
        """
        if not tool_results:
            return "No container information available."

        lines = ["=== CONTAINER & EXECUTION DETAILS (Execute RAG) ===\n"]

        for result in tool_results:
            tool    = result.get("tool", "unknown")
            record  = result.get("record")
            method  = result.get("method", "unknown")

            if record is None:
                lines.append(f"Tool: {tool}\n  ⚠ No container found in BioContainers registry\n")
                continue

            container = record.get("best_container", {})
            uri       = container.get("full_uri", "not available")
            version   = container.get("version", "unknown")

            lines.append(
                f"Tool: {record.get('display_name', tool)}  [lookup: {method}]\n"
                f"  Container URI:  {uri}\n"
                f"  Version:        {version}\n"
                f"  Description:    {record.get('description', 'N/A')[:120]}\n"
                f"  Snakemake:      container: \"{uri}\"\n"
            )

        return "\n".join(lines)
