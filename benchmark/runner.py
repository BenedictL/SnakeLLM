"""
benchmark/runner.py
===================
Orchestrates all benchmark runs across models × RAG conditions × prompts.

Design principles:
  - Cache-first: LLM responses saved as JSON, never re-called if cache exists
  - Resumable: skips runs where output JSON already exists
  - Separates generation (expensive) from evaluation (cheap, re-runnable)
  - Each run produces one JSON in results/runs/ and one row in benchmark.csv

Usage:
    # Run everything
    python -m benchmark.runner

    # Run specific model only
    python -m benchmark.runner --model claude-sonnet-4-6

    # Re-evaluate metrics without re-running LLM calls
    python -m benchmark.runner --eval-only

    # Dry-run: print what would run without calling APIs
    python -m benchmark.runner --dry-run
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

from benchmark.prompts import PROMPTS, BenchmarkPrompt
from benchmark.metrics import evaluate

log = logging.getLogger(__name__)

# ── Configuration ─────────────────────────────────────────────────────────────

MODELS = [
    "claude-sonnet-4-6",
    "gpt-4o",
    "gemini-1.5-pro",
    "qwen2.5-coder-7b",   # fine-tuned
]

RAG_CONDITIONS = [
    "no_rag",
    "plan_rag_only",
    "execute_rag_only",
    "full_dual_rag",
]

RESULTS_DIR  = Path("results")
RUNS_DIR     = RESULTS_DIR / "runs"
CSV_PATH     = RESULTS_DIR / "benchmark.csv"

CSV_FIELDS = [
    # Identity
    "run_id", "model", "rag_condition", "prompt_id",
    "pipeline_type", "difficulty", "timestamp",
    # Metrics
    "valid_spec", "dry_run_pass",
    "rule_count_score", "tool_name_score",
    "container_score", "dag_score", "composite_score",
    # Diagnostics
    "actual_rule_count", "expected_rule_count",
    "actual_tools", "missing_tools",
    "missing_containers", "missing_edges",
    "latency_s", "error",
]


# ── RAG condition builder ─────────────────────────────────────────────────────

def _build_engine(rag_condition: str, model: str):
    """
    Build a SnakeLLMInference engine with the appropriate RAG configuration.
    Returns the engine and a label for logging.
    """
    from llm.plan_rag import PlanRAG
    from llm.execute_rag import ExecuteRAG
    from llm.inference import SnakeLLMInference

    plan_rag    = PlanRAG()    if rag_condition in ("plan_rag_only",    "full_dual_rag") else None
    execute_rag = ExecuteRAG() if rag_condition in ("execute_rag_only", "full_dual_rag") else None

    # Dummy RAG stubs for ablation conditions
    if plan_rag is None:
        plan_rag = _NullRAG()
    if execute_rag is None:
        execute_rag = _NullExecuteRAG()

    return SnakeLLMInference(
        plan_rag=plan_rag,
        execute_rag=execute_rag,
        model=model,
        verbose=False,
    )


class _NullRAG:
    """Stub that returns empty retrieval results — simulates no-RAG condition."""
    def retrieve(self, query: str, top_k: int = 5) -> list:
        return []

    @property
    def collection(self):
        class _C:
            def count(self): return 1  # pretend indexed so no warning
        return _C()


class _NullExecuteRAG:
    """Stub execute RAG — returns empty container URIs."""
    def get_container_uri(self, tool: str) -> str:
        return ""

    def retrieve(self, query: str, top_k: int = 5) -> list:
        return []

    @property
    def collection(self):
        class _C:
            def count(self): return 1
        return _C()


# ── Single run ────────────────────────────────────────────────────────────────

def run_id_for(model: str, rag_condition: str, prompt_id: str) -> str:
    return f"{model}__{rag_condition}__{prompt_id}"


def get_cache_path(run_id: str) -> Path:
    return RUNS_DIR / f"{run_id}.json"


def load_cached_run(run_id: str) -> Optional[dict]:
    path = get_cache_path(run_id)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return None


def save_run(run_id: str, data: dict):
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    get_cache_path(run_id).write_text(
        json.dumps(data, indent=2), encoding="utf-8"
    )


def execute_single_run(
    model:         str,
    rag_condition: str,
    prompt:        BenchmarkPrompt,
    dry_run_flag:  bool = False,
) -> dict:
    """
    Execute one benchmark run. Returns a result dict ready for CSV.

    If cache exists for this run_id, loads it instead of calling the LLM.
    """
    rid = run_id_for(model, rag_condition, prompt.prompt_id)
    cached = load_cached_run(rid)

    if cached and "spec_json" in cached:
        log.info(f"  [CACHE] {rid}")
        spec_json = cached.get("spec_json")
        latency   = cached.get("latency_s", 0)
        error     = cached.get("error", "")
    else:
        if dry_run_flag:
            log.info(f"  [DRY-RUN] would call: {rid}")
            return {"run_id": rid, "model": model, "rag_condition": rag_condition,
                    "prompt_id": prompt.prompt_id, "error": "dry_run_skipped"}

        log.info(f"  [RUN] {rid}")
        spec_json, latency, error = _call_llm(model, rag_condition, prompt)

        save_run(rid, {
            "run_id":        rid,
            "model":         model,
            "rag_condition": rag_condition,
            "prompt_id":     prompt.prompt_id,
            "spec_json":     spec_json,
            "latency_s":     latency,
            "error":         error,
            "timestamp":     datetime.utcnow().isoformat(),
        })

    # Reconstruct PipelineSpec from cached JSON (or None if failed)
    spec = None
    if spec_json:
        try:
            from core.schema import PipelineSpec
            spec = PipelineSpec.model_validate(spec_json)
        except Exception as e:
            error = f"validation_failed: {e}"

    # Evaluate metrics
    metrics = evaluate(spec, prompt, run_dry_run=True)
    metrics.update({
        "run_id":        rid,
        "model":         model,
        "rag_condition": rag_condition,
        "latency_s":     round(latency, 2),
        "timestamp":     datetime.utcnow().isoformat(),
        "error":         error or metrics.get("error", ""),
    })

    return metrics


def _call_llm(
    model: str,
    rag_condition: str,
    prompt: BenchmarkPrompt,
) -> tuple[Optional[dict], float, str]:
    """
    Call the LLM and return (spec_json_dict, latency_seconds, error_string).
    spec_json_dict is None if the call failed.
    """
    try:
        engine = _build_engine(rag_condition, model)
        t0     = time.time()
        spec   = engine.generate(prompt.prompt)
        t1     = time.time()
        return json.loads(spec.model_dump_json()), round(t1 - t0, 2), ""
    except Exception as e:
        return None, 0.0, str(e)


# ── Full benchmark run ────────────────────────────────────────────────────────

def run_benchmark(
    models:         list[str]  = None,
    rag_conditions: list[str]  = None,
    prompts:        list       = None,
    dry_run_flag:   bool       = False,
    eval_only:      bool       = False,
) -> list[dict]:
    """
    Run all combinations, write results to CSV + per-run JSONs.
    Returns list of all result dicts.
    """
    models         = models         or MODELS
    rag_conditions = rag_conditions or RAG_CONDITIONS
    prompts        = prompts        or PROMPTS

    total  = len(models) * len(rag_conditions) * len(prompts)
    done   = 0
    all_results = []

    log.info(f"Benchmark: {len(models)} models × {len(rag_conditions)} RAG conditions × {len(prompts)} prompts = {total} runs")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    RUNS_DIR.mkdir(parents=True, exist_ok=True)

    # Open CSV in append mode so we can resume partial runs
    csv_exists = CSV_PATH.exists()
    csv_file   = open(CSV_PATH, "a", newline="", encoding="utf-8")
    writer     = csv.DictWriter(csv_file, fieldnames=CSV_FIELDS, extrasaction="ignore")
    if not csv_exists:
        writer.writeheader()

    try:
        for model in models:
            for rag_condition in rag_conditions:
                log.info(f"\nModel: {model} | RAG: {rag_condition}")
                for prompt in prompts:
                    rid = run_id_for(model, rag_condition, prompt.prompt_id)

                    # Skip if already in CSV (resume support)
                    if _already_in_csv(rid) and not eval_only:
                        log.info(f"  [SKIP] {rid} already in CSV")
                        done += 1
                        continue

                    if eval_only:
                        # Re-evaluate from cache only — no LLM calls
                        cached = load_cached_run(rid)
                        if not cached:
                            log.warning(f"  [SKIP] no cache for {rid}")
                            continue
                        spec_json = cached.get("spec_json")
                        spec = None
                        if spec_json:
                            try:
                                from core.schema import PipelineSpec
                                spec = PipelineSpec.model_validate(spec_json)
                            except Exception:
                                pass
                        result = evaluate(spec, prompt, run_dry_run=True)
                        result.update({
                            "run_id":        rid,
                            "model":         model,
                            "rag_condition": rag_condition,
                            "latency_s":     cached.get("latency_s", 0),
                            "timestamp":     cached.get("timestamp", ""),
                        })
                    else:
                        result = execute_single_run(model, rag_condition, prompt, dry_run_flag)

                    all_results.append(result)
                    writer.writerow({k: result.get(k, "") for k in CSV_FIELDS})
                    csv_file.flush()

                    done += 1
                    log.info(f"  [{done}/{total}] {rid} → composite={result.get('composite_score', 'N/A'):.3f}")

    finally:
        csv_file.close()

    log.info(f"\nDone. Results → {CSV_PATH}")
    log.info(f"  Per-run JSONs → {RUNS_DIR}/")
    return all_results


def _already_in_csv(run_id: str) -> bool:
    """Check if this run_id already has a row in benchmark.csv."""
    if not CSV_PATH.exists():
        return False
    with open(CSV_PATH, encoding="utf-8") as f:
        return any(run_id in line for line in f)


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(description="SnakeLLM Benchmarking Framework")
    parser.add_argument("--model",      nargs="+", help="Models to run (default: all)")
    parser.add_argument("--rag",        nargs="+", help="RAG conditions (default: all)")
    parser.add_argument("--prompt-id",  nargs="+", help="Specific prompt IDs to run")
    parser.add_argument("--difficulty", choices=["simple", "medium", "hard"], help="Filter by difficulty")
    parser.add_argument("--dry-run",    action="store_true", help="Print runs without calling APIs")
    parser.add_argument("--eval-only",  action="store_true", help="Re-evaluate from cache, no LLM calls")
    args = parser.parse_args()

    # Filter prompts
    selected_prompts = PROMPTS
    if args.prompt_id:
        from benchmark.prompts import PROMPT_BY_ID
        selected_prompts = [PROMPT_BY_ID[pid] for pid in args.prompt_id if pid in PROMPT_BY_ID]
    if args.difficulty:
        selected_prompts = [p for p in selected_prompts if p.difficulty == args.difficulty]

    run_benchmark(
        models         = args.model or MODELS,
        rag_conditions = args.rag   or RAG_CONDITIONS,
        prompts        = selected_prompts,
        dry_run_flag   = args.dry_run,
        eval_only      = args.eval_only,
    )
