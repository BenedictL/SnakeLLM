"""
main.py
=======
SnakeLLM — LLM Core entry point.

Usage:
    # 1. First-time setup (index documents + BioContainers)
    python main.py setup

    # 2. Generate a pipeline from a prompt
    python main.py generate "run differential expression analysis on RNA-seq data using DESeq2"

    # 3. Generate and save to file
    python main.py generate "ATAC-seq peak calling pipeline" --output my_pipeline.json
"""
import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["HF_HUB_DISABLE_PROGRESS_BARS"] = "1"

from dotenv import load_dotenv
load_dotenv()                    # ← this reads your .env file

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import sys
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
if sys.stderr.encoding != 'utf-8':
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

os.makedirs("logs", exist_ok=True)
log_filename = f"logs/snakellm_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.FileHandler(log_filename, encoding="utf-8"),
        logging.StreamHandler()
    ]
)
log = logging.getLogger(__name__)


def check_api_key():
    if os.environ.get("ANTHROPIC_API_KEY"):
        log.info("✓ ANTHROPIC_API_KEY found")
    elif os.environ.get("GEMINI_API_KEY"):
        log.info("✓ GEMINI_API_KEY found")
    elif os.environ.get("OPENAI_API_KEY"):
        log.info("✓ OPENAI_API_KEY found")
    else:
        log.error("No API key found in environment variables.")
        log.error("Please set ANTHROPIC_API_KEY, GEMINI_API_KEY, or OPENAI_API_KEY")
        sys.exit(1)


def cmd_setup(args):
    """
    First-time setup:
    1. Write workflow pattern docs
    2. Fetch BioContainers metadata
    3. Index everything into ChromaDB
    """
    log.info("=" * 60)
    log.info("SnakeLLM — First-Time Setup")
    log.info("=" * 60)

    # Step 1: Generate workflow pattern docs
    log.info("\n[1/3] Writing workflow pattern documents...")
    from data.workflows.workflow_patterns import write_workflow_docs
    write_workflow_docs()

    # Step 2: Fetch BioContainers data
    log.info("\n[2/3] Fetching BioContainers metadata (this may take ~2 min)...")
    from llm.biocontainers_indexer import index_all, DEFAULT_TOOLS
    tools = args.tools if hasattr(args, 'tools') and args.tools else DEFAULT_TOOLS
    results = index_all(tools)
    log.info(f"  ✓ Indexed {len(results)} tools from BioContainers")

    # Step 3: Build RAG indexes
    log.info("\n[3/3] Building RAG vector indexes...")
    from llm.plan_rag import PlanRAG
    from llm.execute_rag import ExecuteRAG

    plan_rag = PlanRAG()
    plan_rag.index_documents(force_reindex=True)

    execute_rag = ExecuteRAG()
    execute_rag.index_documents(force_reindex=True)

    log.info("\n" + "=" * 60)
    log.info("✓ Setup complete! Run: python main.py generate '<your prompt>'")
    log.info("=" * 60)


def cmd_generate(args):
    """Generate a pipeline from a natural language prompt."""
    check_api_key()

    from llm.plan_rag import PlanRAG
    from llm.execute_rag import ExecuteRAG
    from llm.inference import SnakeLLMInference

    log.info(f"\nPrompt: \"{args.prompt}\"")
    log.info("-" * 60)

    # Initialize RAG systems
    plan_rag    = PlanRAG()
    execute_rag = ExecuteRAG()

    # Check if indexed
    if plan_rag.collection.count() == 0 or execute_rag.collection.count() == 0:
        log.warning("RAG indexes are empty. Run: python main.py setup")

    # Run inference
    engine = SnakeLLMInference(
        plan_rag=plan_rag,
        execute_rag=execute_rag,
        verbose=args.verbose if hasattr(args, 'verbose') else False,
    )

    try:
        spec = engine.generate(args.prompt)
    except RuntimeError as e:
        log.error(f"Generation failed: {e}")
        sys.exit(1)

    # Print summary
    log.info("\n" + "=" * 60)
    log.info(f"✓ Pipeline: {spec.pipeline_type}")
    log.info(f"  Description: {spec.description}")
    log.info(f"  Rules ({len(spec.rules)}): {[r.name for r in spec.rules]}")
    log.info(f"  Tools ({len(spec.tools)}): {[t.name for t in spec.tools]}")
    log.info(f"  Execution order: {spec.topological_order()}")
    log.info("=" * 60)

    # Save output
    output_path = args.output if hasattr(args, 'output') and args.output else "pipeline_spec.json"
    spec_json   = json.loads(spec.model_dump_json(indent=2))

    with open(output_path, "w") as f:
        json.dump(spec_json, f, indent=2)

    log.info(f"\n✓ PipelineSpec saved to: {output_path}")
    log.info("  → Hand this file to TM2 (Pipeline Engineer) to generate the Snakefile")

    return spec


def cmd_demo(args):
    """Run a quick demo without API calls to test the RAG pipeline."""
    log.info("Running RAG retrieval demo (no LLM call)...")

    from llm.plan_rag import PlanRAG
    from llm.execute_rag import ExecuteRAG

    plan_rag    = PlanRAG()
    execute_rag = ExecuteRAG()

    test_prompts = [
        "run differential expression analysis on RNA-seq samples",
        "ATAC-seq peak calling with IDR filtering",
        "whole genome sequencing variant calling with GATK",
    ]

    for prompt in test_prompts:
        log.info(f"\nPrompt: '{prompt}'")
        log.info("  Plan RAG hits:")
        plan_hits = plan_rag.retrieve(prompt, top_k=2)
        for h in plan_hits:
            log.info(f"    - {h['metadata'].get('doc_id')} (score: {h.get('score', 0):.3f})")

        # Heuristic tool extraction for demo
        from llm.inference import SnakeLLMInference
        engine = SnakeLLMInference.__new__(SnakeLLMInference)
        engine.plan_rag = plan_rag
        engine.execute_rag = execute_rag
        tools = engine._heuristic_tool_extract(prompt)
        log.info(f"  Detected tools: {tools}")

        log.info("  Execute RAG lookups:")
        for tool in tools[:3]:
            uri = execute_rag.get_container_uri(tool)
            log.info(f"    - {tool}: {uri}")


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SnakeLLM — AI-powered Snakemake pipeline generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command")

    # setup command
    setup_parser = subparsers.add_parser("setup", help="First-time setup: index docs + BioContainers")
    setup_parser.add_argument("--tools", nargs="+", help="Specific tools to index (default: all)")

    # generate command
    gen_parser = subparsers.add_parser("generate", help="Generate a pipeline from a prompt")
    gen_parser.add_argument("prompt",   type=str, help="Natural language pipeline description")
    gen_parser.add_argument("--output", type=str, default="pipeline_spec.json", help="Output JSON path")
    gen_parser.add_argument("--verbose", action="store_true", help="Show raw LLM output")

    # demo command
    demo_parser = subparsers.add_parser("demo", help="Test RAG pipeline without LLM API call")

    args = parser.parse_args()

    if args.command == "setup":
        cmd_setup(args)
    elif args.command == "generate":
        cmd_generate(args)
    elif args.command == "demo":
        cmd_demo(args)
    else:
        parser.print_help()
