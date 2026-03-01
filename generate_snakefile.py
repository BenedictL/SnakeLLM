# generate_snakefile.py
# Usage:  python generate_snakefile.py results/rna_deseq2.json
# Output: output/rna_deseq2/Snakefile + output/rna_deseq2/config.yaml
#
# Changes vs previous version:
#   - clean_shell_cmd() fixes spaces in placeholders + removes bash if/then/fi
#   - Input paths normalised: "data/{sample}_R1.fastq.gz" -> "data/raw/{sample}_R1.fastq.gz"
#   - rule all uses correct expand() for all wildcard outputs
#   - Aggregate rule detection fixed
#   - UTF-8 encoding on all file writes

import json, sys, os, re
from pathlib import Path
from collections import defaultdict, deque

# ── Load pipeline spec ─────────────────────────────────────────────────────────
json_path = sys.argv[1] if len(sys.argv) > 1 else "results/rna_deseq2.json"
spec      = json.load(open(json_path, encoding="utf-8"))

pipeline_name = Path(json_path).stem
base_out      = sys.argv[2] if len(sys.argv) > 2 else "output"
output_dir    = f"{base_out}/{pipeline_name}"
os.makedirs(output_dir, exist_ok=True)

rules     = spec["rules"]
tools_map = {t["name"].lower(): t for t in spec["tools"]}
edges     = spec.get("dag_edges", [])
config    = spec.get("config_params", {})

print(f"Loaded: {pipeline_name}")
print(f"  Rules : {len(rules)}")
print(f"  Tools : {len(spec['tools'])}")
print(f"  Edges : {len(edges)}")

# ── Shell command cleaner ──────────────────────────────────────────────────────
def clean_shell_cmd(cmd):
    """
    Fix all common LLM-generated shell command issues so Snakemake
    can parse the Snakefile without SyntaxError.

    Problems fixed:
      1. Spaced placeholders:  { params . title }  ->  {params.title}
      2. Spaced index:         { input [ 0 ] }     ->  {input[0]}
      3. Spaced simple:        { log }             ->  {log}
      4. Bash if/then/fi conditionals inside $() removed entirely
      5. Duplicate log redirect stripped (we add 2>{log} ourselves)
      6. Double spaces cleaned up
    """

    # ── Step 1: Fix ALL spaced placeholders in one aggressive pass ────────────
    # This handles patterns like { params . title }, { input[0] }, { log[0] }
    # Strategy: find every { ... } block and remove internal spaces

    def fix_placeholder(match):
        inner = match.group(1)
        # Remove all spaces inside
        inner = re.sub(r'\s+', '', inner)
        return '{' + inner + '}'

    # Match { anything } — greedy enough to catch all variants
    # but not so greedy it spans multiple placeholders
    cmd = re.sub(r'\{([^{}]+)\}', fix_placeholder, cmd)

    # ── Step 2: Remove bash conditionals — $( if [ ... ]; then ...; fi ) ──────
    cmd = re.sub(r'\$\(if\s+\[.*?\]\s*;.*?fi\)', '', cmd, flags=re.DOTALL)
    # Also handle: $([ condition ] && echo "flag" || echo "")
    cmd = re.sub(r'\$\(\[.*?\].*?\)', '', cmd, flags=re.DOTALL)

    # ── Step 3: Strip trailing log redirects (we append 2>{log} ourselves) ────
    cmd = re.sub(r'\s*2>\s*\{log\S*\}\s*$', '', cmd).strip()
    cmd = re.sub(r'\s*>\s*\{log\S*\}\s*2>&1\s*$', '', cmd).strip()
    cmd = re.sub(r'\s*&>\s*\{log\S*\}\s*$', '', cmd).strip()

    # ── Step 4: Clean up double spaces left by removals ───────────────────────
    cmd = re.sub(r'  +', ' ', cmd).strip()

    # ── Step 5: Escape ALL quotes so they don't break the Python string literal ─
    cmd = cmd.replace('"', '\\"')

    return cmd

# ── Fix input file paths ───────────────────────────────────────────────────────
def fix_input_path(path):
    """Ensure raw FASTQ inputs always point to data/raw/ directory."""
    # If path starts with data/ but NOT data/raw/, fix it
    if re.match(r'^data/(?!raw/)', path):
        path = path.replace("data/", "data/raw/", 1)
    # If path starts with raw/ without data/, add data/ prefix
    if path.startswith("raw/"):
        path = "data/" + path
    return path

# ── Topological sort ───────────────────────────────────────────────────────────
def topological_sort(rules, edges):
    rule_names = [r["name"] for r in rules]
    in_degree  = {n: 0 for n in rule_names}
    graph      = defaultdict(list)

    for edge in edges:
        src, dst = edge[0], edge[1]
        if src in in_degree and dst in in_degree:
            graph[src].append(dst)
            in_degree[dst] += 1

    queue  = deque([n for n in rule_names if in_degree[n] == 0])
    result = []
    while queue:
        node = queue.popleft()
        result.append(node)
        for neighbor in sorted(graph[node]):
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)

    # Add any rules not reached (disconnected)
    for n in rule_names:
        if n not in result:
            result.append(n)

    return result

ordered_names = topological_sort(rules, edges)
rules_by_name = {r["name"]: r for r in rules}
ordered_rules = [rules_by_name[n] for n in ordered_names if n in rules_by_name]

print(f"  Order : {' -> '.join(ordered_names)}")

# ── Detect terminal rules (never appear as edge source) ───────────────────────
edge_sources      = {e[0] for e in edges}
terminal_names    = [r["name"] for r in rules if r["name"] not in edge_sources]
if not terminal_names:
    terminal_names = [ordered_rules[-1]["name"]]

print(f"  Terminals: {terminal_names}")

# ── Detect aggregate rules (outputs contain no wildcards) ─────────────────────
def is_aggregate(rule):
    return all("{" not in o for o in rule.get("output", []))

def has_wildcard(paths):
    return any("{" in p for p in paths)

# ── Generate rule all targets ──────────────────────────────────────────────────
all_targets = []
for name in terminal_names:
    if name not in rules_by_name:
        continue
    rule = rules_by_name[name]
    for out in rule.get("output", []):
        if "{sample}" in out:
            all_targets.append(f'expand("{out}", sample=config["samples"])')
        else:
            all_targets.append(f'"{out}"')

if not all_targets:
    # Fallback: use last rule in topological order
    last = ordered_rules[-1]
    for out in last.get("output", []):
        if "{sample}" in out:
            all_targets.append(f'expand("{out}", sample=config["samples"])')
        else:
            all_targets.append(f'"{out}"')

# ── Build Snakefile ────────────────────────────────────────────────────────────
lines = []

# Header comment
lines += [
    f'# SnakeLLM — Auto-generated Snakefile',
    f'# Pipeline:    {spec.get("pipeline_type", "unknown")}',
    f'# Description: {spec.get("description", "")}',
    f'# Source:      {json_path}',
    f'',
    'configfile: "config.yaml"',
    '',
    'samples = config["samples"]',
    '',
]

# rule all
lines.append('rule all:')
lines.append('    input:')
for t in all_targets:
    lines.append(f'        {t},')
lines.append('')

# Each rule
for rule in ordered_rules:
    name      = rule["name"]
    inputs    = rule.get("input", [])
    outputs   = rule.get("output", [])
    params    = rule.get("params", {})
    shell_cmd = rule.get("shell_cmd", "")
    script    = rule.get("script", None)
    logs      = rule.get("log", [])
    resources = rule.get("resources", {})
    aggregate = is_aggregate(rule)

    lines.append(f'# {"─" * 60}')
    lines.append(f'rule {name}:')

    # ── input ──
    if inputs:
        lines.append('    input:')
        if aggregate and has_wildcard(inputs):
            # Aggregate rule: wrap wildcard inputs in expand()
            for inp in inputs:
                fixed = fix_input_path(inp)
                if fixed.startswith("expand("):
                    lines.append(f'        {fixed},')
                elif "{sample}" in fixed:
                    lines.append(f'        expand("{fixed}", sample=config["samples"]),')
                else:
                    lines.append(f'        "{fixed}",')
        elif len(inputs) == 1:
            fixed = fix_input_path(inputs[0])
            if fixed.startswith("expand("):
                lines.append(f'        {fixed}')
            else:
                lines.append(f'        "{fixed}"')
        else:
            # Use positional inputs (no labels) to avoid duplicate keyword errors
            # Only use named r1/r2 for exactly 2 paired-end FASTQ inputs
            is_paired_fastq = (
                len(inputs) == 2 and
                all(any(x in inp for x in [".fastq", ".fq", ".gz"]) for inp in inputs)
            )
            if is_paired_fastq:
                for i, inp in enumerate(inputs):
                    fixed = fix_input_path(inp)
                    label = "r1" if i == 0 else "r2"
                    lines.append(f'        {label} = "{fixed}",')
            else:
                # Positional — no labels, avoids duplicate keyword errors entirely
                for inp in inputs:
                    fixed = fix_input_path(inp)
                    if fixed.startswith("expand("):
                        lines.append(f'        {fixed},')
                    else:
                        lines.append(f'        "{fixed}",') 

    # ── output ──
    if outputs:
        lines.append('    output:')
        if len(outputs) == 1:
            lines.append(f'        "{outputs[0]}"')
        else:
            for out in outputs:
                lines.append(f'        "{out}",')

    # ── params ──
    if params:
        lines.append('    params:')
        for k, v in params.items():
            if isinstance(v, str):
                lines.append(f'        {k} = "{v}",')
            else:
                lines.append(f'        {k} = {v},')

    # ── threads ──
    cpus = resources.get("cpus", 1)
    lines.append(f'    threads: {cpus}')

    # ── resources ──
    mem  = resources.get("mem_mb", 4000)
    time = resources.get("time_min", 60)
    lines.append('    resources:')
    lines.append(f'        mem_mb  = {mem},')
    lines.append(f'        runtime = {time}')

    # ── log ──
    if logs:
        lines.append('    log:')
        lines.append(f'        "{logs[0]}"')
    else:
        log_path = f'logs/{name}/{{sample}}.log' if not aggregate else f'logs/{name}/{name}.log'
        lines.append('    log:')
        lines.append(f'        "{log_path}"')

    # ── shell or script ──
    if script and not shell_cmd:
        lines.append('    script:')
        lines.append(f'        "{script}"')
    elif shell_cmd:
        cleaned = clean_shell_cmd(shell_cmd)
        lines.append('    shell:')
        lines.append(f'        "{cleaned} 2> {{log}}"')
    else:
        lines.append('    shell:')
        lines.append(f'        "echo \'{name} — no shell command defined\'"')

    lines.append('')

snakefile_content = "\n".join(lines)

# ── Generate config.yaml ───────────────────────────────────────────────────────

# Extract sample list from config_params FIRST, then remove it to avoid duplication
sample_list = config.get("samples", ["sample1", "sample2", "sample3"])
if not isinstance(sample_list, list):
    sample_list = ["sample1", "sample2", "sample3"]
config_no_samples = {k: v for k, v in config.items() if k != "samples"}

config_lines = [
    f'# SnakeLLM — Pipeline configuration',
    f'# Pipeline: {spec.get("pipeline_type", "unknown")}',
    f'# Edit paths and parameters before running',
    '',
    '# Sample names — edit to match your actual sample IDs',
    'samples:',
]
for s in sample_list:
    config_lines.append(f'  - {s}')
config_lines.append('')

# Categorise remaining config params (samples already removed)
ref_keys   = [k for k in config_no_samples if any(x in k.lower() for x in
              ["dir","file","path","gtf","fa","fasta","index","adapter","genome","annot","ref"])]
stat_keys  = [k for k in config_no_samples if any(x in k.lower() for x in
              ["threshold","cutoff","alpha","pval","padj","lfc","min","max","filter"])]
other_keys = [k for k in config_no_samples if k not in ref_keys and k not in stat_keys]

if ref_keys:
    config_lines.append('# Reference files — update paths for your system')
    for k in ref_keys:
        v = config_no_samples[k]
        config_lines.append(f'{k}: "{v}"' if isinstance(v, str) else f'{k}: {v}')
    config_lines.append('')

if stat_keys:
    config_lines.append('# Statistical thresholds')
    for k in stat_keys:
        v = config_no_samples[k]
        config_lines.append(f'{k}: {v}')
    config_lines.append('')

if other_keys:
    config_lines.append('# Pipeline parameters')
    for k in other_keys:
        v = config_no_samples[k]
        config_lines.append(f'{k}: "{v}"' if isinstance(v, str) else f'{k}: {v}')
    config_lines.append('')

config_content = "\n".join(config_lines)

# ── Write files ────────────────────────────────────────────────────────────────
snakefile_path = f"{output_dir}/Snakefile"
config_path    = f"{output_dir}/config.yaml"

with open(snakefile_path, "w", encoding="utf-8") as f:
    f.write(snakefile_content)

with open(config_path, "w", encoding="utf-8") as f:
    f.write(config_content)

# ── Create dummy input files for dry-run ──────────────────────────────────────
raw_dir = f"{output_dir}/data/raw"
os.makedirs(raw_dir, exist_ok=True)

# Use the actual sample list from the pipeline spec
for s in sample_list:
    for suffix in ["_R1.fastq.gz", "_R2.fastq.gz"]:
        dummy = f"{raw_dir}/{s}{suffix}"
        if not os.path.exists(dummy):
            open(dummy, "w", encoding="utf-8").close()

# Create common QC directories needed for MultiQC dry-run
qc_dirs = [
    "logs/star", "logs/bowtie2", "logs/fastqc", "logs/trimmomatic",
    "logs/hisat2", "logs/bwa", "logs/picard",
    "qc/raw", "qc/trimmed", "qc/raw_fastqc", "qc/post_trim_fastqc",
    "qc/fastqc", "qc/flagstat", "aligned", "counts", "results",
]
for d in qc_dirs:
    os.makedirs(f"{output_dir}/{d}", exist_ok=True)

print(f"\nGenerated:")
print(f"  {snakefile_path}")
print(f"  {config_path}")
print(f"  {raw_dir}/sample[1-3]_R[12].fastq.gz  (dummy files for dry-run)")
print(f"\nTo dry-run:")
print(f"  cd {output_dir}")
print(f"  snakemake --dry-run --cores 1")