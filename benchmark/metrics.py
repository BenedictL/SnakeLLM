"""
benchmark/metrics.py
====================
Metric evaluation for a single benchmark run.

Given a PipelineSpec (or None if generation failed) and a BenchmarkPrompt,
computes all 6 metrics and returns a flat dict ready for benchmark.csv.

Metrics:
  1. valid_spec        (0/1) — did Pydantic validation pass?
  2. dry_run_pass      (0/1) — does snakemake --dry-run exit 0?
  3. rule_count_score  (0.0–1.0) — actual / expected rule count, capped at 1.0
  4. tool_name_score   (0.0–1.0) — fraction of expected tools present
  5. container_score   (0.0–1.0) — fraction of expected container substrings present
  6. dag_score         (0.0–1.0) — fraction of expected edges present in topo order

All scores are floats in [0, 1]. The composite score is their unweighted mean.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from benchmark.prompts import BenchmarkPrompt


def evaluate(
    spec,                        # PipelineSpec | None
    prompt: BenchmarkPrompt,
    run_dry_run: bool = True,    # set False in unit tests or if snakemake not installed
) -> dict:
    """
    Evaluate a generation result against the benchmark ground truth.

    Args:
        spec:        Validated PipelineSpec from inference.py, or None if failed
        prompt:      The BenchmarkPrompt used to generate this spec
        run_dry_run: Whether to invoke snakemake --dry-run (requires snakemake)

    Returns:
        Flat dict of metric scores, ready to append to benchmark.csv
    """
    result = {
        "prompt_id":        prompt.prompt_id,
        "pipeline_type":    prompt.pipeline_type,
        "difficulty":       prompt.difficulty,
        # scores
        "valid_spec":       0,
        "dry_run_pass":     0,
        "rule_count_score": 0.0,
        "tool_name_score":  0.0,
        "container_score":  0.0,
        "dag_score":        0.0,
        "composite_score":  0.0,
        # diagnostics
        "actual_rule_count":    0,
        "expected_rule_count":  prompt.expected_rule_count,
        "actual_tools":         "",
        "missing_tools":        "",
        "missing_containers":   "",
        "missing_edges":        "",
        "error":                "",
    }

    if spec is None:
        result["error"] = "generation_failed"
        return result

    # ── 1. valid_spec ─────────────────────────────────────────────────────────
    result["valid_spec"] = 1  # if we got here, Pydantic already validated it

    # ── 2. rule_count_score ───────────────────────────────────────────────────
    actual_count = len(spec.rules)
    result["actual_rule_count"] = actual_count
    result["rule_count_score"]  = min(actual_count / prompt.expected_rule_count, 1.0)

    # ── 3. tool_name_score ────────────────────────────────────────────────────
    actual_tools   = {t.name.lower() for t in spec.tools}
    expected_tools = {t.lower() for t in prompt.expected_tools}
    # Fuzzy match: expected tool is "present" if any actual tool name contains it
    # (handles bwa-mem2 vs bwamem2 vs bwa_mem2 variants)
    matched_tools  = {
        exp for exp in expected_tools
        if any(_tool_matches(exp, act) for act in actual_tools)
    }
    missing_tools  = expected_tools - matched_tools
    result["actual_tools"]      = ", ".join(sorted(actual_tools))
    result["missing_tools"]     = ", ".join(sorted(missing_tools))
    result["tool_name_score"]   = len(matched_tools) / len(expected_tools) if expected_tools else 1.0

    # ── 4. container_score ────────────────────────────────────────────────────
    # Collect all container URIs from tools
    container_uris = []
    for tool in spec.tools:
        if hasattr(tool.container, "full_uri"):
            container_uris.append(tool.container.full_uri.lower())
        elif isinstance(tool.container, str):
            container_uris.append(tool.container.lower())

    all_uris = " ".join(container_uris)
    matched_containers = [
        sub for sub in prompt.expected_containers
        if sub.lower() in all_uris
    ]
    missing_containers = [
        sub for sub in prompt.expected_containers
        if sub.lower() not in all_uris
    ]
    result["missing_containers"] = ", ".join(missing_containers)
    result["container_score"]    = (
        len(matched_containers) / len(prompt.expected_containers)
        if prompt.expected_containers else 1.0
    )

    # ── 5. dag_score ──────────────────────────────────────────────────────────
    if prompt.expected_dag_edges:
        actual_edges   = {(src.lower(), dst.lower()) for src, dst in spec.dag_edges}
        rule_name_map  = {r.name.lower(): r.name.lower() for r in spec.rules}

        matched_edges  = []
        missing_edges  = []

        for exp_src, exp_dst in prompt.expected_dag_edges:
            # Fuzzy: check if any actual edge contains the expected src/dst as substrings
            found = any(
                exp_src in act_src and exp_dst in act_dst
                for act_src, act_dst in actual_edges
            )
            if found:
                matched_edges.append((exp_src, exp_dst))
            else:
                missing_edges.append(f"{exp_src}→{exp_dst}")

        result["missing_edges"] = ", ".join(missing_edges)
        result["dag_score"]     = len(matched_edges) / len(prompt.expected_dag_edges)
    else:
        result["dag_score"] = 1.0  # no edges expected (single-rule pipeline)

    # ── 6. dry_run_pass ───────────────────────────────────────────────────────
    if run_dry_run and shutil.which("snakemake"):
        result["dry_run_pass"] = _run_dry_run(spec)
    else:
        result["dry_run_pass"] = -1  # -1 = not evaluated

    # ── composite score ───────────────────────────────────────────────────────
    scores = [
        result["valid_spec"],
        result["rule_count_score"],
        result["tool_name_score"],
        result["container_score"],
        result["dag_score"],
    ]
    # Only include dry_run in composite if it was evaluated
    if result["dry_run_pass"] != -1:
        scores.append(result["dry_run_pass"])

    result["composite_score"] = sum(scores) / len(scores)

    return result


def _tool_matches(expected: str, actual: str) -> bool:
    """
    Fuzzy tool name matching to handle registry naming inconsistencies.

    Examples:
      bwa-mem2  matches  bwamem2, bwa_mem2, bwa-mem2
      trim-galore matches trimgalore, trim_galore
    """
    # Normalise both to alphanumeric only
    norm_exp = "".join(c for c in expected.lower() if c.isalnum())
    norm_act = "".join(c for c in actual.lower() if c.isalnum())
    return norm_exp in norm_act or norm_act in norm_exp


def _run_dry_run(spec) -> int:
    """
    Assemble the spec into a temp dir and run snakemake --dry-run.
    Returns 1 if exit code 0, else 0.
    """
    try:
        from pipeline.assembler import assemble
        with tempfile.TemporaryDirectory() as tmpdir:
            snakefile, _ = assemble(
                spec,
                output_dir=tmpdir,
                source_name="benchmark",
                create_dummy_inputs=True,
            )
            result = subprocess.run(
                ["snakemake", "--dry-run", "--cores", "1",
                 "--snakefile", str(snakefile),
                 "--directory", tmpdir],
                capture_output=True,
                text=True,
                timeout=60,
            )
            return 1 if result.returncode == 0 else 0
    except Exception:
        return 0
