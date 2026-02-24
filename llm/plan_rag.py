"""
llm/plan_rag.py
===============
PLAN RAG — Stage 1 retrieval
"What tools should I use, in what order, with what dependencies?"

Retrieves high-level workflow patterns and tool relationship knowledge.
Uses hybrid retrieval: semantic (ChromaDB) + keyword (BM25) + RRF fusion.

This answers the PLANNING question before the LLM generates the PipelineSpec.
"""

from __future__ import annotations
import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
from rank_bm25 import BM25Okapi
import chromadb
from chromadb.utils import embedding_functions
from sentence_transformers import SentenceTransformer

log = logging.getLogger(__name__)

WORKFLOW_DOCS_DIR  = Path("data/workflows")
TOOL_DOCS_DIR      = Path("data/tools")
CHROMA_DIR         = Path("data/chroma/plan_rag")
EMBED_MODEL        = "allenai-specter"   # good for scientific text; fallback: all-MiniLM-L6-v2


class PlanRAG:
    """
    Retrieves workflow patterns and tool-relationship knowledge to answer:
    "Given this analysis request, what is the correct pipeline structure?"

    Knowledge base contains:
    - Workflow patterns (rna_seq_de.md, atac_seq.md, wgs.md, ...)
    - Tool relationship docs  (which tools pair well, common DAG patterns)
    - File format chain docs  (FASTQ → BAM → VCF etc.)
    """

    def __init__(self, embed_model: str = EMBED_MODEL):
        log.info(f"Initializing Plan RAG (embed model: {embed_model})")
        self.embed_model_name = embed_model
        self._load_embedding_model()
        self._load_chromadb()
        self._load_bm25_index()

    # ── INITIALIZATION ────────────────────────────────────────────────────────

    def _load_embedding_model(self):
        """Load sentence transformer for dense embeddings."""
        try:
            self.embed_model = SentenceTransformer(self.embed_model_name)
            log.info(f"  ✓ Loaded embedding model: {self.embed_model_name}")
        except Exception as e:
            log.warning(f"  Could not load {self.embed_model_name}, falling back to all-MiniLM-L6-v2: {e}")
            self.embed_model = SentenceTransformer("all-MiniLM-L6-v2")

    def _load_chromadb(self):
        """Load or create the ChromaDB vector store for Plan RAG."""
        CHROMA_DIR.mkdir(parents=True, exist_ok=True)
        self.chroma_client = chromadb.PersistentClient(path=str(CHROMA_DIR))
        self.collection    = self.chroma_client.get_or_create_collection(
            name="plan_rag",
            metadata={"hnsw:space": "cosine"}
        )
        count = self.collection.count()
        log.info(f"  ✓ Plan RAG ChromaDB loaded ({count} documents indexed)")

    def _load_bm25_index(self):
        """Build in-memory BM25 index from all plan documents."""
        self.bm25_docs     = []
        self.bm25_metadata = []

        for source_dir in [WORKFLOW_DOCS_DIR, TOOL_DOCS_DIR]:
            if not source_dir.exists():
                continue
            for f in source_dir.glob("*.md"):
                self.bm25_docs.append(f.read_text())
                self.bm25_metadata.append({
                    "doc_id":   f.stem,
                    "source":   source_dir.name,
                    "filepath": str(f)
                })

        if self.bm25_docs:
            tokenized    = [doc.lower().split() for doc in self.bm25_docs]
            self.bm25    = BM25Okapi(tokenized)
            log.info(f"  ✓ BM25 index built ({len(self.bm25_docs)} documents)")
        else:
            self.bm25 = None
            log.warning("  No documents found for BM25 index. Run index_documents() first.")

    # ── INDEXING ──────────────────────────────────────────────────────────────

    def index_documents(self, force_reindex: bool = False):
        """
        (Re)index all workflow and tool documents into ChromaDB.
        Call this once after adding or updating docs.
        """
        if self.collection.count() > 0 and not force_reindex:
            log.info("Plan RAG already indexed. Pass force_reindex=True to rebuild.")
            return

        log.info("Indexing Plan RAG documents...")
        docs, ids, metadatas = [], [], []

        for source_dir in [WORKFLOW_DOCS_DIR, TOOL_DOCS_DIR]:
            if not source_dir.exists():
                log.warning(f"  Directory not found: {source_dir}")
                continue
            for f in source_dir.glob("*.md"):
                text = f.read_text().strip()
                docs.append(text)
                ids.append(f"{source_dir.name}::{f.stem}")
                metadatas.append({
                    "doc_id": f.stem,
                    "source": source_dir.name,
                    "category": self._extract_category(text),
                })

        if not docs:
            log.error("No documents to index.")
            return

        # Embed in batches
        batch_size = 32
        all_embeddings = []
        for i in range(0, len(docs), batch_size):
            batch = docs[i:i+batch_size]
            embeddings = self.embed_model.encode(batch, show_progress_bar=False).tolist()
            all_embeddings.extend(embeddings)

        # Upsert into ChromaDB
        self.collection.upsert(
            documents=docs,
            embeddings=all_embeddings,
            ids=ids,
            metadatas=metadatas,
        )
        log.info(f"  ✓ Indexed {len(docs)} documents into Plan RAG")

        # Rebuild BM25
        self._load_bm25_index()

    def _extract_category(self, text: str) -> str:
        """Extract category tag from document header."""
        for line in text.split("\n"):
            if line.startswith("Category:"):
                return line.replace("Category:", "").strip()
        return "unknown"

    # ── RETRIEVAL ─────────────────────────────────────────────────────────────

    def retrieve(self, prompt: str, top_k: int = 4) -> list[dict]:
        """
        Main retrieval entry point.
        Returns top-k most relevant workflow/tool documents via hybrid RRF search.
        """
        expanded = self._expand_query(prompt)

        semantic_hits = self._semantic_search(expanded, top_k=10)
        bm25_hits     = self._bm25_search(expanded, top_k=10)
        fused         = self._rrf_fusion(semantic_hits, bm25_hits, final_top_k=top_k)

        log.debug(f"Plan RAG: {len(fused)} docs retrieved for prompt: '{prompt[:60]}...'")
        return fused

    def _semantic_search(self, query: str, top_k: int) -> list[dict]:
        """Dense vector search via ChromaDB."""
        if self.collection.count() == 0:
            log.warning("Plan RAG collection is empty. Run index_documents() first.")
            return []

        embedding = self.embed_model.encode(query).tolist()
        results   = self.collection.query(
            query_embeddings=[embedding],
            n_results=min(top_k, self.collection.count()),
            include=["documents", "metadatas", "distances"]
        )

        return [
            {
                "doc":      results["documents"][0][i],
                "metadata": results["metadatas"][0][i],
                "score":    1 - results["distances"][0][i],  # distance → similarity
                "method":   "semantic"
            }
            for i in range(len(results["documents"][0]))
        ]

    def _bm25_search(self, query: str, top_k: int) -> list[dict]:
        """Sparse BM25 keyword search."""
        if self.bm25 is None or not self.bm25_docs:
            return []

        tokens = query.lower().split()
        scores = self.bm25.get_scores(tokens)
        top_idx = np.argsort(scores)[::-1][:top_k]

        return [
            {
                "doc":      self.bm25_docs[i],
                "metadata": self.bm25_metadata[i],
                "score":    float(scores[i]),
                "method":   "bm25"
            }
            for i in top_idx if scores[i] > 0
        ]

    def _rrf_fusion(
        self,
        semantic_hits: list[dict],
        bm25_hits: list[dict],
        k: int = 60,
        final_top_k: int = 4
    ) -> list[dict]:
        """
        Reciprocal Rank Fusion — merges semantic and BM25 rankings.
        Score = Σ 1 / (k + rank_in_list)
        k=60 is the standard constant that prevents top-ranked docs from dominating.
        """
        scores    = {}   # doc_id → fused RRF score
        doc_store = {}   # doc_id → full doc dict

        for rank, hit in enumerate(semantic_hits):
            doc_id = hit["metadata"].get("doc_id", f"sem_{rank}")
            scores[doc_id]    = scores.get(doc_id, 0) + 1 / (k + rank + 1)
            doc_store[doc_id] = hit

        for rank, hit in enumerate(bm25_hits):
            doc_id = hit["metadata"].get("doc_id", f"bm25_{rank}")
            scores[doc_id]    = scores.get(doc_id, 0) + 1 / (k + rank + 1)
            doc_store[doc_id] = hit

        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        return [doc_store[doc_id] for doc_id, _ in ranked[:final_top_k]]

    def _expand_query(self, prompt: str) -> str:
        """
        Synonym expansion to improve recall for bioinformatics terminology.
        Adds related terms that might appear in the knowledge base docs.
        """
        synonyms = {
            "differential expression":  "DEG DESeq2 edgeR DE analysis transcriptomics gene expression",
            "de analysis":               "differential expression DESeq2 edgeR RNA-seq",
            "rna-seq":                   "RNA sequencing transcriptomics alignment quantification",
            "rna seq":                   "RNA sequencing transcriptomics alignment quantification",
            "atac-seq":                  "chromatin accessibility open chromatin MACS2 peak calling Tn5",
            "atac seq":                  "chromatin accessibility open chromatin MACS2 peak calling",
            "chip-seq":                  "histone modification peak calling MACS2 ChIP immunoprecipitation",
            "wgs":                       "whole genome sequencing variant calling GATK SNP indel",
            "variant calling":           "SNP indel GATK HaplotypeCaller VCF genotyping WGS",
            "peak calling":              "MACS2 ATAC-seq ChIP-seq open chromatin accessible",
            "single cell":               "scRNA-seq 10x Genomics Cell Ranger Seurat Scanpy clustering",
            "methylation":               "WGBS bisulfite Bismark CpG methylation epigenomics",
            "alignment":                 "mapping STAR BWA HISAT2 read alignment",
            "trimming":                  "adapter removal Trimmomatic fastp quality trimming",
            "quality control":           "QC FastQC MultiQC trimming adapter",
            "enrichment":                "GO enrichment KEGG pathway clusterProfiler functional annotation",
        }

        expanded = [prompt]
        prompt_lower = prompt.lower()
        for key, expansion in synonyms.items():
            if key in prompt_lower:
                expanded.append(expansion)

        return " ".join(expanded)

    # ── CONTEXT FORMATTING ────────────────────────────────────────────────────

    def format_for_prompt(self, retrieved_docs: list[dict]) -> str:
        """
        Formats retrieved Plan RAG docs for injection into the LLM system prompt.
        Returns clean, structured text that the LLM can reason over.
        """
        if not retrieved_docs:
            return "No relevant workflow patterns found."

        sections = ["=== WORKFLOW & TOOL KNOWLEDGE (Plan RAG) ===\n"]
        for i, hit in enumerate(retrieved_docs, 1):
            source  = hit["metadata"].get("source", "unknown")
            doc_id  = hit["metadata"].get("doc_id", f"doc_{i}")
            score   = hit.get("score", 0)
            sections.append(
                f"--- Document {i}: {doc_id} (source: {source}, relevance: {score:.3f}) ---\n"
                f"{hit['doc']}\n"
            )

        return "\n".join(sections)
