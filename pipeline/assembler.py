"""
pipeline/assembler.py
=====================
Converts a validated PipelineSpec into a runnable Snakefile + config.yaml.

Owned by: Pipeline Engineer (TM2 — Rinoshan)
Consumes: core.schema.PipelineSpec  (produced by LLM Core — Benedict)

This is the integration point between Layer 2 (LLM Core) and Layer 3
(Pipeline Generation). The PipelineSpec is already validated by Pydantic
before this module receives it — no dict manipulation, no regex patching.

Usage:
    from core.schema import PipelineSpec
    from pipeline.assembler import assemble, assemble_from_file

    # From a PipelineSpec object (normal path — output of inference.py)
    spec = ...  # already validated PipelineSpec
    assemble(spec, output_dir="output/rna_deseq2")

    # From a saved JSON file (dev/testing path)
    assemble_from_file("results/rna_deseq2.json", output_dir="output/rna_deseq2")
"""

from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader

from core.schema import PipelineSpec, RuleSpec

log = logging.getLogger(__name__)

# Templates live alongside the assembler in pipeline/templates/
# Falls back to root-level templates/ for backward compatibility
_TEMPLATE_DIRS = [
    Path(__file__).parent / "templates",   # pipeline/templates/ (correct location)
    Path(__file__).parent.parent / "templates",  # root templates/ (backward compat)
]


# ── TEMPLATE ENVIRONMENT ──────────────────────────────────────────────────────

def _get_jinja_env() -> Environment:
    """Load Jinja2 environment from the first template directory that exists."""
    for d in _TEMPLATE_DIRS:
        if d.exists():
            log.debug(f"Loading templates from: {d}")
            return Environment(loader=FileSystemLoader(str(d)))
    raise FileNotFoundError(
        f"No templates directory found. Checked: {[str(d) for d in _TEMPLATE_DIRS]}"
    )


# ── SHELL COMMAND SANITISER ───────────────────────────────────────────────────

def _sanitise_shell_cmd(cmd: str) -> str:
    """
    Light sanitisation of LLM-generated shell commands.

    Only fixes structural issues that would cause Snakemake SyntaxError:
    - Collapses whitespace inside {wildcard} placeholders
    - Strips trailing log redirects (the template adds 2> {log} itself)
    - Collapses double spaces

    Does NOT attempt to fix wrong paths or wrong flags — those are prompt
    engineering problems, not assembler problems.
    """
    if not cmd:
        return cmd

    # Collapse whitespace inside { } placeholders: { sample } → {sample}
    cmd = re.sub(r'\{([^{}]+)\}', lambda m: '{' + re.sub(r'\s+', '', m.group(1)) + '}', cmd)

    # Strip trailing log redirects — snakefile.j2 adds "2> {log}" itself
    cmd = re.sub(r'\s*2>\s*\{log[^}]*\}\s*$', '', cmd).strip()
    cmd = re.sub(r'\s*>\s*\{log[^}]*\}\s*2>&1\s*$', '', cmd).strip()
    cmd = re.sub(r'\s*&>\s*\{log[^}]*\}\s*$', '', cmd).strip()

    # Collapse double spaces
    cmd = re.sub(r'  +', ' ', cmd).strip()

    # Escape internal double quotes for safe embedding in Snakemake shell: "..."
    cmd = cmd.replace('"', '\\"')

    return cmd


# ── INPUT PATH NORMALISER ─────────────────────────────────────────────────────

def _normalise_input_path(path: str) -> str:
    """
    Ensure raw FASTQ inputs always reference data/raw/.

    This corrects the two most common LLM path hallucinations:
      data/sample.fastq.gz  →  data/raw/sample.fastq.gz
      raw/sample.fastq.gz   →  data/raw/sample.fastq.gz

    All other paths are returned unchanged.
    """
    if re.match(r'^data/(?!raw/)', path):
        return path.replace("data/", "data/raw/", 1)
    if path.startswith("raw/"):
        return "data/" + path
    return path


# ── RULE TEMPLATE DATA BUILDER ────────────────────────────────────────────────

def _is_aggregate_rule(rule: RuleSpec) -> bool:
    """
    An aggregate rule produces outputs with no wildcards — e.g. MultiQC,
    DESeq2. These need expand() on their inputs instead of wildcard patterns.
    """
    return all("{" not in o for o in rule.output)


def _has_wildcard(paths: list[str]) -> bool:
    return any("{" in p for p in paths)


def _build_rule_template_data(rule: RuleSpec) -> dict[str, Any]:
    """
    Converts a validated RuleSpec into a flat dict for Jinja2 rendering.
    Uses spec fields directly — no dict.get() fallbacks needed.
    """
    aggregate = _is_aggregate_rule(rule)

    # ── Process input paths ───────────────────────────────────────────────────
    processed_inputs: list[str] = []

    if aggregate and _has_wildcard(rule.input):
        # Aggregate rule with wildcard inputs → wrap in expand()
        for inp in rule.input:
            fixed = _normalise_input_path(inp)
            if fixed.startswith("expand("):
                processed_inputs.append(fixed)
            elif "{sample}" in fixed:
                processed_inputs.append(f'expand("{fixed}", sample=config["samples"])')
            else:
                processed_inputs.append(f'"{fixed}"')

    elif len(rule.input) == 2 and all("fastq" in p for p in rule.input):
        # Paired-end FASTQ — name as r1/r2 (template handles this)
        for inp in rule.input:
            processed_inputs.append(f'"{_normalise_input_path(inp)}"')

    else:
        for inp in rule.input:
            fixed = _normalise_input_path(inp)
            if fixed.startswith("expand("):
                processed_inputs.append(fixed)
            else:
                processed_inputs.append(f'"{fixed}"')

    # ── Log path ──────────────────────────────────────────────────────────────
    if rule.log:
        log_path = rule.log[0]
    elif aggregate:
        log_path = f"logs/{rule.name}/{rule.name}.log"
    else:
        log_path = f"logs/{rule.name}/{{sample}}.log"

    # ── Threads: use RuleSpec.threads if set, else fall back to resources.cpus ─
    threads = rule.threads if rule.threads is not None else rule.resources.cpus

    return {
        "name":      rule.name,
        "inputs":    processed_inputs,
        "outputs":   rule.output,
        "params":    rule.params,
        "shell_cmd": _sanitise_shell_cmd(rule.shell_cmd) if rule.shell_cmd else "",
        "script":    rule.script,
        "log":       log_path,
        "threads":   threads,
        "resources": {
            "mem_mb":   rule.resources.mem_mb,
            "time_min": rule.resources.time_min,
            "disk_mb":  rule.resources.disk_mb,
        },
    }


# ── ALL-RULE TARGETS ─────────────────────────────────────────────────────────

def _build_all_targets(spec: PipelineSpec) -> list[str]:
    """
    Determine the final outputs for rule all: by finding terminal rules
    (rules that are not sources of any DAG edge).
    Uses spec.topological_order() — no reimplementation.
    """
    edge_sources   = {src for src, _ in spec.dag_edges}
    terminal_names = [r.name for r in spec.rules if r.name not in edge_sources]

    # Fallback: last rule in topological order
    if not terminal_names:
        order          = spec.topological_order()
        terminal_names = [order[-1]] if order else []

    rules_by_name = {r.name: r for r in spec.rules}
    targets: list[str] = []

    for name in terminal_names:
        rule = rules_by_name.get(name)
        if not rule:
            continue
        for out in rule.output:
            if "{sample}" in out:
                targets.append(f'expand("{out}", sample=config["samples"])')
            else:
                targets.append(f'"{out}"')

    return targets


# ── CONFIG PARAM CATEGORISER ──────────────────────────────────────────────────

_REF_KEYWORDS  = {"dir", "file", "path", "gtf", "fa", "fasta", "index",
                  "adapter", "genome", "annot", "ref", "blacklist", "bed", "dict"}
_STAT_KEYWORDS = {"threshold", "cutoff", "alpha", "pval", "padj", "lfc",
                  "min", "max", "filter", "fdr", "fc"}


def _categorise_config(config: dict[str, Any]) -> tuple[dict, dict, dict]:
    """
    Split config_params into three groups for cleaner config.yaml output:
    - ref_keys:   file paths, genome references
    - stat_keys:  statistical thresholds
    - other_keys: everything else (pipeline params, booleans, sample names)
    """
    no_samples = {k: v for k, v in config.items() if k != "samples"}

    ref_keys   = {k: v for k, v in no_samples.items()
                  if any(kw in k.lower() for kw in _REF_KEYWORDS)}
    stat_keys  = {k: v for k, v in no_samples.items()
                  if k not in ref_keys and any(kw in k.lower() for kw in _STAT_KEYWORDS)}
    other_keys = {k: v for k, v in no_samples.items()
                  if k not in ref_keys and k not in stat_keys}

    return ref_keys, stat_keys, other_keys


# ── CORE ASSEMBLY ─────────────────────────────────────────────────────────────

def render(spec: PipelineSpec, source_name: str = "unknown") -> tuple[str, str]:
    """
    Render a PipelineSpec into (snakefile_content, config_content) strings.

    This is the pure rendering function — no file I/O, no side effects.
    Useful for testing and for the API layer (TM4).

    Args:
        spec:        Validated PipelineSpec object
        source_name: Label shown in Snakefile header comment

    Returns:
        (snakefile_str, config_yaml_str)
    """
    env = _get_jinja_env()

    # Build ordered rule list using PipelineSpec's own topological sort
    order         = spec.topological_order()
    rules_by_name = {r.name: r for r in spec.rules}
    ordered_rules = [rules_by_name[n] for n in order if n in rules_by_name]

    # Build template context
    samples                      = spec.config_params.get("samples", ["sample1", "sample2", "sample3"])
    if not isinstance(samples, list):
        samples = ["sample1", "sample2", "sample3"]

    ref_keys, stat_keys, other_keys = _categorise_config(spec.config_params)

    context = {
        "spec":        spec,
        "source_name": source_name,
        "all_targets": _build_all_targets(spec),
        "rules":       [_build_rule_template_data(r) for r in ordered_rules],
        "samples":     samples,
        "ref_keys":    ref_keys,
        "stat_keys":   stat_keys,
        "other_keys":  other_keys,
    }

    snakefile_content = env.get_template("snakefile.j2").render(**context)
    config_content    = env.get_template("config.j2").render(**context)

    return snakefile_content, config_content


def assemble(
    spec:       PipelineSpec,
    output_dir: str,
    source_name: str = "unknown",
    create_dummy_inputs: bool = True,
) -> tuple[Path, Path]:
    """
    Full assembly: render PipelineSpec → write Snakefile + config.yaml to disk.

    Args:
        spec:                Validated PipelineSpec
        output_dir:          Directory to write Snakefile and config.yaml into
        source_name:         Label for Snakefile header comment
        create_dummy_inputs: If True, create empty dummy FASTQ files for dry-run

    Returns:
        (snakefile_path, config_path)
    """
    snakefile_content, config_content = render(spec, source_name)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    snakefile_path = out / "Snakefile"
    config_path    = out / "config.yaml"

    snakefile_path.write_text(snakefile_content, encoding="utf-8")
    config_path.write_text(config_content, encoding="utf-8")

    log.info(f"  ✓ Snakefile written: {snakefile_path}")
    log.info(f"  ✓ config.yaml written: {config_path}")

    if create_dummy_inputs:
        _create_dummy_inputs(spec, out)

    log.info(f"\nTo dry-run:")
    log.info(f"  cd {output_dir}")
    log.info(f"  snakemake --dry-run --cores 1")

    return snakefile_path, config_path


def assemble_from_file(
    json_path:   str,
    output_dir:  str = None,
    base_out:    str = "output",
) -> tuple[Path, Path]:
    """
    Convenience wrapper: load a PipelineSpec JSON → validate → assemble.

    This is the dev/testing path — equivalent to the old generate_snakefile.py
    but now goes through proper PipelineSpec validation first.

    Args:
        json_path:  Path to a saved PipelineSpec JSON file
        output_dir: Explicit output directory (overrides base_out)
        base_out:   Base output directory (output_dir = base_out/pipeline_name)
    """
    path = Path(json_path)
    if not path.exists():
        raise FileNotFoundError(f"PipelineSpec JSON not found: {json_path}")

    raw = json.loads(path.read_text(encoding="utf-8"))

    log.info(f"Loading: {path.name}")
    spec = PipelineSpec.model_validate(raw)  # validates on load — no silent failures
    log.info(f"  Rules : {len(spec.rules)}")
    log.info(f"  Tools : {len(spec.tools)}")
    log.info(f"  Edges : {len(spec.dag_edges)}")

    out_dir = output_dir or f"{base_out}/{path.stem}"
    return assemble(spec, out_dir, source_name=str(json_path))


# ── DUMMY INPUT CREATOR (for dry-run testing) ─────────────────────────────────

_DUMMY_LOG_DIRS = [
    "logs/star", "logs/bowtie2", "logs/fastqc", "logs/trimmomatic",
    "logs/hisat2", "logs/bwa", "logs/picard", "logs/featurecounts",
    "logs/umi_dedup", "logs/deseq2", "logs/macs2", "logs/samtools",
    "qc/raw", "qc/trimmed", "qc/raw_fastqc", "qc/post_trim_fastqc",
    "qc/fastqc", "qc/flagstat", "aligned", "counts", "results",
]


def _create_dummy_inputs(spec: PipelineSpec, out: Path):
    """
    Create empty dummy FASTQ files and log directories so snakemake --dry-run
    doesn't fail on missing input files.
    """
    samples = spec.config_params.get("samples", ["sample1", "sample2", "sample3"])
    if not isinstance(samples, list):
        samples = ["sample1", "sample2", "sample3"]

    raw_dir = out / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        for suffix in ["_R1.fastq.gz", "_R2.fastq.gz"]:
            dummy = raw_dir / f"{sample}{suffix}"
            if not dummy.exists():
                dummy.touch()

    for d in _DUMMY_LOG_DIRS:
        (out / d).mkdir(parents=True, exist_ok=True)

    log.info(f"  ✓ Dummy inputs created in {raw_dir}")


# ── CLI (backward-compatible with generate_snakefile.py) ──────────────────────

if __name__ == "__main__":
    import sys
    import logging as _logging
    _logging.basicConfig(level=_logging.INFO, format="%(levelname)s  %(message)s")

    json_path  = sys.argv[1] if len(sys.argv) > 1 else "results/rna_deseq2.json"
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None

    snakefile, config = assemble_from_file(json_path, output_dir=output_dir)
    print(f"\nGenerated:")
    print(f"  {snakefile}")
    print(f"  {config}")