# generate_snakefile.py
# Usage:  python generate_snakefile.py results/rna_deseq2.json
# Output: output/rna_deseq2/Snakefile + output/rna_deseq2/config.yaml

import json
import sys
import os
import re
from pathlib import Path
from collections import defaultdict, deque
from jinja2 import Environment, FileSystemLoader

# ── Shell command cleaner ──────────────────────────────────────────────────────
def clean_shell_cmd(cmd):
    """
    Fix all common LLM-generated shell command issues so Snakemake
    can parse the Snakefile without SyntaxError.
    """
    def fix_placeholder(match):
        inner = match.group(1)
        inner = re.sub(r'\s+', '', inner)
        return '{' + inner + '}'

    cmd = re.sub(r'\{([^{}]+)\}', fix_placeholder, cmd)
    cmd = re.sub(r'\$\(if\s+\[.*?\]\s*;.*?fi\)', '', cmd, flags=re.DOTALL)
    cmd = re.sub(r'\$\(\[.*?\].*?\)', '', cmd, flags=re.DOTALL)
    cmd = re.sub(r'\s*2>\s*\{log\S*\}\s*$', '', cmd).strip()
    cmd = re.sub(r'\s*>\s*\{log\S*\}\s*2>&1\s*$', '', cmd).strip()
    cmd = re.sub(r'\s*&>\s*\{log\S*\}\s*$', '', cmd).strip()
    cmd = re.sub(r'  +', ' ', cmd).strip()
    cmd = cmd.replace('"', '\\"')
    return cmd

# ── Fix input file paths ───────────────────────────────────────────────────────
def fix_input_path(path):
    """Ensure raw FASTQ inputs always point to data/raw/ directory."""
    if re.match(r'^data/(?!raw/)', path):
        path = path.replace("data/", "data/raw/", 1)
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

    for n in rule_names:
        if n not in result:
            result.append(n)

    return result

def is_aggregate(rule):
    return all("{" not in o for o in rule.get("output", []))

def has_wildcard(paths):
    return any("{" in p for p in paths)

# ── Generator functions ────────────────────────────────────────────────────────

def prepare_template_data(spec: dict, source_name: str) -> dict:
    """Prepare all the derived variables needed by the Jinja2 templates."""
    rules     = spec.get("rules", [])
    edges     = spec.get("dag_edges", [])
    config    = spec.get("config_params", {})

    ordered_names = topological_sort(rules, edges)
    rules_by_name = {r["name"]: r for r in rules}
    ordered_rules = [rules_by_name[n] for n in ordered_names if n in rules_by_name]

    edge_sources      = {e[0] for e in edges}
    terminal_names    = [r["name"] for r in rules if r["name"] not in edge_sources]
    if not terminal_names and ordered_rules:
        terminal_names = [ordered_rules[-1]["name"]]

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

    if not all_targets and ordered_rules:
        last = ordered_rules[-1]
        for out in last.get("output", []):
            if "{sample}" in out:
                all_targets.append(f'expand("{out}", sample=config["samples"])')
            else:
                all_targets.append(f'"{out}"')

    # Pre-process rule fields for the template
    template_rules = []
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

        processed_inputs = []
        if aggregate and has_wildcard(inputs):
            for inp in inputs:
                fixed = fix_input_path(inp)
                if fixed.startswith("expand("):
                    processed_inputs.append(fixed)
                elif "{sample}" in fixed:
                    processed_inputs.append(f'expand("{fixed}", sample=config["samples"])')
                else:
                    processed_inputs.append(f'"{fixed}"')
        elif len(inputs) == 1:
            fixed = fix_input_path(inputs[0])
            if fixed.startswith("expand("):
                processed_inputs.append(fixed)
            else:
                processed_inputs.append(f'"{fixed}"')
        else:
            for inp in inputs:
                fixed = fix_input_path(inp)
                if fixed.startswith("expand("):
                    processed_inputs.append(fixed)
                else:
                    processed_inputs.append(f'"{fixed}"')

        log_path = logs[0] if logs else f'logs/{name}/{{sample}}.log' if not aggregate else f'logs/{name}/{name}.log'
        cleaned_shell = clean_shell_cmd(shell_cmd) if shell_cmd else ""

        template_rules.append({
            "name": name,
            "inputs": processed_inputs,
            "outputs": outputs,
            "params": params,
            "resources": resources,
            "log": log_path,
            "script": script,
            "shell_cmd": cleaned_shell
        })

    # Prepare config.yaml variables
    sample_list = config.get("samples", ["sample1", "sample2", "sample3"])
    if not isinstance(sample_list, list):
        sample_list = ["sample1", "sample2", "sample3"]
    config_no_samples = {k: v for k, v in config.items() if k != "samples"}

    ref_keys   = {k: v for k, v in config_no_samples.items() if any(x in k.lower() for x in
                  ["dir","file","path","gtf","fa","fasta","index","adapter","genome","annot","ref"])}
    stat_keys  = {k: v for k, v in config_no_samples.items() if any(x in k.lower() for x in
                  ["threshold","cutoff","alpha","pval","padj","lfc","min","max","filter"])}
    other_keys = {k: v for k, v in config_no_samples.items() if k not in ref_keys and k not in stat_keys}

    return {
        "spec": spec,
        "source_name": source_name,
        "all_targets": all_targets,
        "rules": template_rules,
        "samples": sample_list,
        "ref_keys": ref_keys,
        "stat_keys": stat_keys,
        "other_keys": other_keys
    }

def generate_contents(spec: dict, source_name: str = "unknown") -> tuple[str, str]:
    """Generates the content for Snakefile and config.yaml using Jinja2 templates."""
    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), "templates")))
    
    template_data = prepare_template_data(spec, source_name)
    
    snakefile_template = env.get_template("snakefile.j2")
    config_template = env.get_template("config.j2")
    
    snakefile_content = snakefile_template.render(**template_data)
    config_content = config_template.render(**template_data)
    
    return snakefile_content, config_content


def generate_pipeline(spec: dict, output_dir: str, source_name: str = "unknown"):
    """Generates the pipeline files and creates dry-run dummy directories."""
    snakefile_content, config_content = generate_contents(spec, source_name)

    os.makedirs(output_dir, exist_ok=True)
    snakefile_path = os.path.join(output_dir, "Snakefile")
    config_path    = os.path.join(output_dir, "config.yaml")

    with open(snakefile_path, "w", encoding="utf-8") as f:
        f.write(snakefile_content)
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(config_content)

    # Create dummy input files for dry-run
    raw_dir = os.path.join(output_dir, "data", "raw")
    os.makedirs(raw_dir, exist_ok=True)

    sample_list = spec.get("config_params", {}).get("samples", ["sample1", "sample2", "sample3"])
    if not isinstance(sample_list, list):
        sample_list = ["sample1", "sample2", "sample3"]

    for s in sample_list:
        for suffix in ["_R1.fastq.gz", "_R2.fastq.gz"]:
            dummy = os.path.join(raw_dir, f"{s}{suffix}")
            if not os.path.exists(dummy):
                open(dummy, "w", encoding="utf-8").close()

    qc_dirs = [
        "logs/star", "logs/bowtie2", "logs/fastqc", "logs/trimmomatic",
        "logs/hisat2", "logs/bwa", "logs/picard", "logs/featurecounts",
        "logs/umi_dedup", "logs/deseq2",
        "qc/raw", "qc/trimmed", "qc/raw_fastqc", "qc/post_trim_fastqc",
        "qc/fastqc", "qc/flagstat", "aligned", "counts", "results",
    ]
    for d in qc_dirs:
        os.makedirs(os.path.join(output_dir, d), exist_ok=True)

    print(f"\nGenerated:")
    print(f"  {snakefile_path}")
    print(f"  {config_path}")
    print(f"  {raw_dir}/sample[1-3]_R[12].fastq.gz  (dummy files for dry-run)")
    print(f"\nTo dry-run:")
    print(f"  cd {output_dir}")
    print(f"  snakemake --dry-run --cores 1")


def main():
    json_path = sys.argv[1] if len(sys.argv) > 1 else "results/rna_deseq2.json"
    if not os.path.exists(json_path):
        print(f"Error: Could not find {json_path}")
        sys.exit(1)
        
    with open(json_path, encoding="utf-8") as f:
        spec = json.load(f)

    pipeline_name = Path(json_path).stem
    base_out      = sys.argv[2] if len(sys.argv) > 2 else "output"
    output_dir    = f"{base_out}/{pipeline_name}"

    print(f"Loaded: {pipeline_name}")
    print(f"  Rules : {len(spec.get('rules', []))}")
    print(f"  Tools : {len(spec.get('tools', []))}")
    print(f"  Edges : {len(spec.get('dag_edges', []))}")

    generate_pipeline(spec, output_dir, source_name=json_path)

if __name__ == "__main__":
    main()