"""
benchmark/analyse.py
====================
Summarise benchmark.csv results into publication-ready tables.

Usage:
    python -m benchmark.analyse                    # full summary
    python -m benchmark.analyse --by difficulty    # breakdown by difficulty
    python -m benchmark.analyse --model claude-sonnet-4-6  # single model deep-dive
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

CSV_PATH = Path("results/benchmark.csv")

METRICS = [
    "valid_spec", "dry_run_pass", "rule_count_score",
    "tool_name_score", "container_score", "dag_score", "composite_score",
]


def load_results() -> list[dict]:
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"No results found at {CSV_PATH}. Run benchmark first.")
    with open(CSV_PATH, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def summarise_by(results: list[dict], group_keys: list[str]) -> dict:
    """Group results and compute mean metrics per group."""
    groups = defaultdict(list)
    for row in results:
        key = tuple(row[k] for k in group_keys)
        groups[key].append(row)

    summary = {}
    for key, rows in sorted(groups.items()):
        group_label = " | ".join(f"{k}={v}" for k, v in zip(group_keys, key))
        summary[group_label] = {
            "n": len(rows),
            **{
                metric: round(_mean([
                    float(r[metric]) for r in rows
                    if r.get(metric) not in ("", "-1", None)
                    and r.get(metric) != "-1"
                ]), 3)
                for metric in METRICS
            }
        }
    return summary


def print_table(summary: dict, title: str):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

    # Header
    col_w = 14
    header = f"{'Group':<35} {'N':>4}"
    for m in METRICS:
        short = m.replace("_score", "").replace("_", " ")[:col_w]
        header += f"  {short:>{col_w}}"
    print(header)
    print("-" * len(header))

    for group, metrics in summary.items():
        row = f"{group:<35} {metrics['n']:>4}"
        for m in METRICS:
            val = metrics.get(m, 0.0)
            row += f"  {val:>{col_w}.3f}"
        print(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--by",     default="model", choices=["model", "rag_condition", "difficulty", "pipeline_type"])
    parser.add_argument("--model",  help="Filter to single model")
    parser.add_argument("--rag",    help="Filter to single RAG condition")
    args = parser.parse_args()

    results = load_results()
    print(f"\nLoaded {len(results)} benchmark runs from {CSV_PATH}")

    # Apply filters
    if args.model:
        results = [r for r in results if r["model"] == args.model]
    if args.rag:
        results = [r for r in results if r["rag_condition"] == args.rag]

    # Main summary table
    summary = summarise_by(results, [args.by])
    print_table(summary, f"Results by {args.by}")

    # Always show model × RAG condition matrix
    if args.by != "model":
        model_rag = summarise_by(results, ["model", "rag_condition"])
        print_table(model_rag, "Model × RAG condition matrix")

    # Difficulty breakdown
    diff_summary = summarise_by(results, ["difficulty"])
    print_table(diff_summary, "Results by difficulty")

    # Composite score ranking
    print(f"\n{'='*80}")
    print("  Composite Score Ranking (model × RAG)")
    print(f"{'='*80}")
    model_rag_full = summarise_by(results, ["model", "rag_condition"])
    ranked = sorted(model_rag_full.items(), key=lambda x: x[1]["composite_score"], reverse=True)
    for i, (group, metrics) in enumerate(ranked, 1):
        print(f"  {i:2}. {group:<50} composite={metrics['composite_score']:.3f}  n={metrics['n']}")


if __name__ == "__main__":
    main()
