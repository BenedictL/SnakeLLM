"""
llm/inference.py
================
The LLM inference layer — ties Plan RAG + Execute RAG into one generation call.

Flow:
  1. Plan RAG   → retrieve workflow patterns + tool relationships
  2. LLM Call 1 → generate tool plan (which tools, in what order)
  3. Execute RAG → retrieve exact container URIs for planned tools
  4. LLM Call 2 → generate full PipelineSpec JSON (grounded with containers)
  5. Validate   → pydantic schema check + retry loop on failure

This two-call architecture is intentional:
- Call 1 is cheap (small output, planning only)
- Call 2 has full grounding from Execute RAG
- Separating them gives us a natural retry boundary
"""

from __future__ import annotations
import json
import logging
import re
from typing import Optional
from pydantic import ValidationError
import json_repair

from core.schema import PipelineSpec
from llm.plan_rag import PlanRAG
from llm.execute_rag import ExecuteRAG
from llm.providers import get_provider, LLMProvider

log = logging.getLogger(__name__)

# ── PROMPTS ───────────────────────────────────────────────────────────────────

PLAN_SYSTEM_PROMPT = """You are a senior bioinformatics pipeline architect.
Given a natural language analysis request, identify:
1. The analysis type (e.g. rna-seq-de, atac-seq, wgs-variant-calling)
2. The ordered list of tools needed
3. The pipeline dependencies

Use the workflow knowledge provided. Think step by step.

Respond with ONLY this JSON (no markdown, no prose):
{{
  "pipeline_type": "<type>",
  "description": "<one sentence>",
  "ordered_tools": ["tool1", "tool2", ...],
  "tool_purposes": {{"tool1": "purpose", "tool2": "purpose"}},
  "dag_edges": [["tool1_rule", "tool2_rule"], ...]
}}

WORKFLOW KNOWLEDGE:
{plan_context}
"""

EXECUTE_SYSTEM_PROMPT = """You are a Snakemake pipeline code generator.
Generate a complete PipelineSpec JSON for the given analysis.

RULES:
1. Output ONLY valid JSON matching the exact schema below. No markdown, no prose, no explanation.
2. Use the EXACT container URIs from the Container Details section — do NOT invent URIs.
3. shell_cmd must be valid shell syntax with Snakemake {{wildcards}} and {{params}} substitutions.
4. dag_edges must reference actual rule names you define.
5. Think about input/output file format chaining before writing rules.
6. IMPORTANT: Rules that produce a single output file (no wildcards in output) must use aggregate inputs. DESeq2, MultiQC, and clusterProfiler are aggregate rules — their input must reference a merged file or directory, NOT a {{sample}} wildcard pattern.
7. strandedness must always be a config_param, never hardcoded. Include in config_params:
  "strandedness": 0,  # 0=unstranded, 1=stranded, 2=reverse-stranded
  "paired_end": true  # true for PE, false for SE
The featureCounts -p flag must only be used when paired_end is true.

PIPELINESPEC SCHEMA:
{schema}

ANALYSIS PLAN:
{plan}

CONTAINER DETAILS (use these EXACT URIs):
{execute_context}

FEW-SHOT EXAMPLES:
{examples}
"""

# ── FEW-SHOT EXAMPLES ─────────────────────────────────────────────────────────
# These teach the LLM the exact output format. Add more for better accuracy.

FEW_SHOT_EXAMPLES = """
EXAMPLE 1 — Single rule (trim_reads):
Input prompt: "trim adapter sequences from paired-end RNA-seq reads"
Output:
{
  "pipeline_type": "preprocessing",
  "description": "Adapter trimming of paired-end RNA-seq reads using Trimmomatic",
  "tools": [
    {
      "name": "Trimmomatic",
      "version": "0.39",
      "container": {
        "registry": "quay.io",
        "image": "biocontainers/trimmomatic",
        "tag": "0.39--hdfd78af_2",
        "full_uri": "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2",
        "source": "biocontainers"
      },
      "purpose": "Remove adapter sequences and low-quality bases from reads",
      "language": "CLI"
    }
  ],
  "rules": [
    {
      "name": "trim_reads",
      "tool": "Trimmomatic",
      "input": ["data/{{sample}}_R1.fastq.gz", "data/{{sample}}_R2.fastq.gz"],
      "output": ["trimmed/{{sample}}_R1.fastq.gz", "trimmed/{{sample}}_R2.fastq.gz"],
      "params": {"adapters": "{{config[adapters]}}", "min_len": 36, "threads": 4},
      "shell_cmd": "trimmomatic PE -threads {{params.threads}} {{input[0]}} {{input[1]}} {{output[0]}} /dev/null {{output[1]}} /dev/null ILLUMINACLIP:{{params.adapters}}:2:30:10 MINLEN:{{params.min_len}}",
      "resources": {"cpus": 4, "mem_mb": 8000, "time_min": 60, "disk_mb": 20000},
      "log": ["logs/trimmomatic/{{sample}}.log"]
    },
    {
      "name": "align_star",
      "tool": "STAR",
      "input": ["trimmed/{{sample}}_R1.fastq.gz", "trimmed/{{sample}}_R2.fastq.gz"],
      "output": ["aligned/{{sample}}.bam"],
      "params": {"genomeDir": "{{config[genome_build]}}"},
      "shell_cmd": "STAR --runThreadN 8 --genomeDir {{params.genomeDir}} --readFilesIn {{input[0]}} {{input[1]}} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix aligned/{{sample}}_",
      "resources": {"cpus": 8, "mem_mb": 32000, "time_min": 120, "disk_mb": 50000},
      "log": ["logs/star/{{sample}}.log"]
    },
    {
      "name": "run_deseq2",
      "tool": "DESeq2",
      "input": ["counts/all_samples_counts.txt"],
      "output": ["results/diff_expr.csv"],
      "params": {},
      "shell_cmd": "Rscript scripts/run_deseq2.R {{input[0]}} {{output[0]}}",
      "resources": {"cpus": 4, "mem_mb": 8000, "time_min": 60, "disk_mb": 20000},
      "log": ["logs/deseq2/diff_expr.log"]
    }
  ],
  "dag_edges": [
    ["trim_reads", "post_trim_fastqc"],
    ["trim_reads", "align_star"],
    ["align_star", "run_deseq2"]
  ],
  "config_params": {
    "adapters": "resources/adapters.fa", 
    "samples": ["S1", "S2"],
    "strandedness": 2,
    "paired_end": true,
    "organism_db": "org.Hs.eg.db",
    "genome_build": "hg38",
    "gtf_source": "ensembl"
  },
  "wildcards": ["sample"]
}
"""


# ── MAIN INFERENCE CLASS ──────────────────────────────────────────────────────

class SnakeLLMInference:
    """
    Orchestrates the full Dual RAG + two-call LLM generation pipeline.
    """

    def __init__(
        self,
        plan_rag:    PlanRAG,
        execute_rag: ExecuteRAG,
        provider:    LLMProvider = None,
        max_retries: int = 3,
        verbose:     bool = False,
    ):
        self.plan_rag    = plan_rag
        self.execute_rag = execute_rag
        self.provider    = provider or get_provider()
        self.max_retries = max_retries
        self.verbose     = verbose

        if verbose:
            logging.basicConfig(level=logging.DEBUG)

    # ── MAIN ENTRY POINT ─────────────────────────────────────────────────────

    def generate(self, user_prompt: str) -> PipelineSpec:
        """
        Full generation pipeline:
        1. Plan RAG retrieval
        2. LLM planning call (what tools?)
        3. Execute RAG retrieval (exact containers for planned tools)
        4. LLM generation call (full PipelineSpec)
        5. Pydantic validation + retry loop
        """
        log.info(f"Generating pipeline for: '{user_prompt}' using {self.provider.name}")

        # ── STEP 1: Plan RAG ─────────────────────────────────────────────────
        log.info("  [1/4] Plan RAG retrieval...")
        plan_docs = self.plan_rag.retrieve(user_prompt, top_k=4)
        plan_context = self.plan_rag.format_for_prompt(plan_docs)

        # ── STEP 2: LLM Planning Call ────────────────────────────────────────
        log.info("  [2/4] LLM planning call...")
        tool_plan = self._planning_call(user_prompt, plan_context)
        log.info(f"         → Tools identified: {tool_plan.get('ordered_tools', [])}")

        # ── STEP 3: Execute RAG ──────────────────────────────────────────────
        log.info("  [3/4] Execute RAG container lookup...")
        tool_names     = tool_plan.get("ordered_tools", [])
        execute_hits   = self.execute_rag.retrieve_for_tools(tool_names)
        execute_context = self.execute_rag.format_for_prompt(execute_hits)

        # ── STEP 4: LLM Generation Call with Retry ───────────────────────────
        log.info("  [4/4] LLM generation call (PipelineSpec)...")
        spec = self._generation_call_with_retry(
            user_prompt=user_prompt,
            tool_plan=tool_plan,
            execute_context=execute_context,
        )

        log.info(f"  ✓ Pipeline generated: {len(spec.rules)} rules, {len(spec.tools)} tools")
        return spec

    # ── LLM CALLS ─────────────────────────────────────────────────────────────

    def _planning_call(self, user_prompt: str, plan_context: str) -> dict:
        """
        Call 1 — lightweight planning call.
        Extracts: pipeline type, ordered tools, DAG edges.
        """
        system = PLAN_SYSTEM_PROMPT.format(plan_context=plan_context)

        try:
            raw = self.provider.complete(
                system=system,
                messages=[{"role": "user", "content": user_prompt}],
                max_tokens=1024,
            ).strip()
            raw  = self._strip_markdown_fences(raw)
            plan = json_repair.loads(raw)
            return plan
        except Exception as e:
            log.warning(f"Planning call failed ({e}), using fallback plan")
            # Fallback: extract tool names from prompt heuristically
            return {
                "pipeline_type": "unknown",
                "description":   user_prompt,
                "ordered_tools": self._heuristic_tool_extract(user_prompt),
                "dag_edges":     []
            }

    def _generation_call_with_retry(
        self,
        user_prompt: str,
        tool_plan: dict,
        execute_context: str,
    ) -> PipelineSpec:
        """
        Call 2 — full PipelineSpec generation with validation retry loop.
        On schema failure, feeds the exact pydantic error back to the LLM.
        """
        schema_json = PipelineSpec.model_json_schema()
        system = EXECUTE_SYSTEM_PROMPT.format(
            schema=json.dumps(schema_json, indent=2),
            plan=json.dumps(tool_plan, indent=2),
            execute_context=execute_context,
            examples=FEW_SHOT_EXAMPLES,
        )

        messages = [{"role": "user", "content": user_prompt}]
        last_good_json = None  # tracks last JSON that parsed but failed schema

        for attempt in range(1, self.max_retries + 1):
            log.info(f"    Generation attempt {attempt}/{self.max_retries}...")

            raw_text = self.provider.complete(
                system=system,
                messages=[messages[0]],  # always send only the original user message
                max_tokens=8192,
            ).strip()

            if self.verbose:
                log.debug(f"Raw LLM output:\n{raw_text[:500]}...")

            # Try to parse JSON
            try:
                raw_json = json_repair.loads(self._strip_markdown_fences(raw_text))
                last_good_json = raw_json  # save the last parseable JSON
            except ValueError as e:
                log.warning(f"    JSON parse failed: {e}")
                continue

            # Try to validate schema
            try:
                spec = PipelineSpec(**raw_json)
                return spec  # ✓ success
            except ValidationError as e:
                errors = self._format_validation_errors(e)
                log.warning(f"    Schema validation failed:\n{errors}")
                print("EXACT PYDANTIC ERROR:", e.errors())  # ← add this

                # Retry with focused correction prompt (don't grow the history)
                system += (
                    f"\n\nPREVIOUS ATTEMPT ERRORS (fix these):\n{errors}"
                )

        # All attempts exhausted — save the last parseable JSON
        import os
        os.makedirs("results", exist_ok=True)
        raw_path = "results/my_pipeline.json"
        if last_good_json:
            with open(raw_path, "w", encoding="utf-8") as f:
                json.dump(last_good_json, f, indent=2)
            log.info(f"  ✓ Raw pipeline saved to: {raw_path} (schema validation skipped)")
            raise RuntimeError(
                f"Validation failed after {self.max_retries} attempts. "
                f"Raw output saved to {raw_path}."
            )
        raise RuntimeError(
            f"LLM failed to produce valid JSON after {self.max_retries} attempts."
        )

    # ── HELPERS ───────────────────────────────────────────────────────────────

    def _strip_markdown_fences(self, text: str) -> str:
        """Remove ```json ... ``` fences from LLM output."""
        text = re.sub(r"^```(?:json)?\s*", "", text.strip(), flags=re.MULTILINE)
        text = re.sub(r"```\s*$", "", text.strip(), flags=re.MULTILINE)
        return text.strip()

    def _format_validation_errors(self, e: ValidationError) -> str:
        """Format pydantic ValidationError into a clear message for the LLM."""
        lines = []
        for err in e.errors():
            loc  = " → ".join(str(x) for x in err["loc"])
            msg  = err["msg"]
            lines.append(f"  • {loc}: {msg}")
        return "\n".join(lines)

    def _heuristic_tool_extract(self, prompt: str) -> list[str]:
        """
        Last-resort: extract known tool names mentioned in the prompt.
        Used only when the planning call fails entirely.
        """
        known_tools = [
            "star", "hisat2", "bwa", "bwa-mem2", "bowtie2",
            "trimmomatic", "fastp", "fastqc", "multiqc",
            "featurecounts", "htseq", "samtools", "picard",
            "deseq2", "edger", "clusterprofiler",
            "macs2", "macs3", "deeptools",
            "gatk4", "bcftools", "deepvariant",
            "bismark", "cellranger",
        ]
        prompt_lower = prompt.lower()
        return [t for t in known_tools if t in prompt_lower]