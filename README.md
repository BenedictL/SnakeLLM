# SnakeLLM — LLM Core

> AI-powered Snakemake pipeline generator — LLM Core module  

---

## Architecture: Dual RAG

```
User Prompt
    │
    ├──────────────────────────────────────┐
    ▼                                      ▼
PLAN RAG                             EXECUTE RAG
"What tools, in what order?"        "Exact container URI + CLI?"
  ↓ Workflow patterns                  ↓ BioContainers registry
  ↓ Tool relationships                 ↓ Pinned image tags
  ↓ File format chains                 ↓ Resource benchmarks
    │                                      │
    └──────────────┬───────────────────────┘
                   ▼
            LLM (Claude)
                   ▼
           PipelineSpec JSON
           (validated by pydantic)
                   ▼
            → TM2 (Pipeline Engineer)
              generates Snakefile
```

---

## Quick Start

### 1. Install dependencies
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Set your API key
```bash
export ANTHROPIC_API_KEY=sk-ant-...
```

### 3. First-time setup (run once)
```bash
# Generates workflow docs + fetches BioContainers + builds vector indexes
python main.py setup
```

### 4. Generate a pipeline
```bash
python main.py generate "run differential expression analysis on paired-end RNA-seq data using DESeq2"

python main.py generate "ATAC-seq chromatin accessibility pipeline with peak calling and motif analysis"

python main.py generate "WGS variant calling with GATK4 best practices for 10 samples"
```

### 5. Test RAG without API call
```bash
python main.py demo
```

---

## Project Structure

```
snakellm/
├── core/
│   └── schema.py              # PipelineSpec pydantic models (shared with TM2)
├── llm/
│   ├── plan_rag.py            # Plan RAG: workflow patterns + tool relationships
│   ├── execute_rag.py         # Execute RAG: BioContainers container URIs
│   ├── biocontainers_indexer.py  # Fetches + indexes BioContainers registry
│   └── inference.py           # Main generation engine (two-call LLM pipeline)
├── data/
│   ├── workflows/             # Workflow pattern .md docs (Plan RAG knowledge base)
│   ├── tools/                 # Tool profile .md docs (Plan RAG)
│   ├── biocontainers/         # BioContainers JSON records (Execute RAG)
│   └── chroma/                # ChromaDB vector stores (auto-generated)
├── tests/
├── main.py                    # CLI entry point
└── requirements.txt
```

---

## Output

`generate` produces a `pipeline_spec.json` — the PipelineSpec schema agreed with TM2:

```json
{
  "pipeline_type": "rna-seq-de",
  "description": "...",
  "tools": [...],
  "rules": [...],
  "dag_edges": [["trim_reads", "align_star"], ...],
  "config_params": {...},
  "wildcards": ["sample"]
}
```

Hand this file to TM2. They consume it to generate the Snakefile.

---

## Re-indexing BioContainers

To add new tools to the Execute RAG knowledge base:
```bash
# Add specific tools
python -m llm.biocontainers_indexer --tools cellranger spaceranger

# Re-index from scratch
python main.py setup
```

---

## Evaluation

To measure RAG quality (schema pass rate, tool accuracy):
```bash
python tests/benchmark.py
```
