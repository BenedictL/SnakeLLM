"""
core/schema.py
==============
The PipelineSpec is the central contract between the LLM Core (Benedict)
and the Pipeline Engineer (TM2). Every field here is agreed upon by both.

NEVER change field names without notifying TM2 — their Jinja2 templates
depend on this schema directly.
"""

from __future__ import annotations
from typing import Any, Optional
from pydantic import BaseModel, Field, field_validator


# ── TOOL LEVEL ────────────────────────────────────────────────────────────────

class ContainerRef(BaseModel):
    """A pinned, reproducible container image for a single tool."""
    registry:   str = Field(..., description="e.g. quay.io, docker.io, ghcr.io")
    image:      str = Field(..., description="e.g. biocontainers/star")
    tag:        str = Field(..., description="e.g. 2.7.10a--h9ee0642_1")
    full_uri:   str = Field(..., description="full pullable URI: registry/image:tag")
    source:     str = Field(default="biocontainers", description="biocontainers | custom | bioconductor")


class ToolSpec(BaseModel):
    """One tool used in the pipeline."""
    name:        str           = Field(..., description="canonical tool name e.g. STAR")
    version:     str           = Field(..., description="pinned version e.g. 2.7.10a")
    container:   ContainerRef
    purpose:     str           = Field(..., description="one-line description of role in pipeline")
    language:    str           = Field(default="CLI", description="CLI | R | Python")


# ── RULE LEVEL ────────────────────────────────────────────────────────────────

class ResourceSpec(BaseModel):
    """Computational resource requirements for a Snakemake rule."""
    cpus:       int   = Field(default=4,     ge=1)
    mem_mb:     int   = Field(default=8000,  ge=256)
    time_min:   int   = Field(default=120,   ge=1,   description="wall time in minutes")
    disk_mb:    int   = Field(default=10000, ge=256)


class RuleSpec(BaseModel):
    """
    One Snakemake rule. Maps 1-to-1 with a rule block in the generated Snakefile.
    TM2 uses this to fill Jinja2 templates.
    """
    name:        str             = Field(..., description="snake_case rule name")
    tool:        str             = Field(..., description="must match a ToolSpec.name")
    input:       list[str]       = Field(..., description="input file patterns with {wildcards}")
    output:      list[str]       = Field(..., description="output file patterns with {wildcards}")
    params:      dict[str, Any]  = Field(default_factory=dict, description="key-value CLI params")
    shell_cmd:   str             = Field(..., description="shell command template")
    script:      Optional[str]   = Field(default=None, description="path to R/Python script if not shell")
    resources:   ResourceSpec    = Field(default_factory=ResourceSpec)
    log:         list[str]       = Field(default_factory=list, description="log file paths")

    @field_validator("name")
    @classmethod
    def name_is_snake_case(cls, v: str) -> str:
        import re
        if not re.match(r'^[a-z][a-z0-9_]*$', v):
            raise ValueError(f"Rule name '{v}' must be snake_case")
        return v


# ── PIPELINE LEVEL ────────────────────────────────────────────────────────────

class PipelineSpec(BaseModel):
    """
    The complete structured pipeline specification.
    Produced by:  LLM Core (Benedict)
    Consumed by:  Pipeline Engineer (TM2) → generates Snakefile + config.yaml
                  Infrastructure (TM3)    → adds container directives
    """
    pipeline_type:  str                  = Field(..., description="e.g. rna-seq-de, atac-seq, wgs-variant")
    description:    str                  = Field(..., description="one sentence summary of what this pipeline does")
    tools:          list[ToolSpec]        = Field(..., min_length=1)
    rules:          list[RuleSpec]        = Field(..., min_length=1)
    dag_edges:      list[tuple[str, str]] = Field(..., description="(rule_from, rule_to) dependency pairs")
    config_params:  dict[str, Any]        = Field(..., description="user-tunable config.yaml defaults")
    wildcards:      list[str]             = Field(default_factory=list, description="e.g. ['sample', 'unit']")

    @field_validator("dag_edges")
    @classmethod
    def edges_reference_valid_rules(cls, edges: list, info) -> list:
        # Only validate if rules are already parsed
        if hasattr(info, 'data') and 'rules' in info.data:
            rule_names = {r.name.lower() for r in info.data["rules"]}
            for src, dst in edges:
                if src.lower() not in rule_names:
                    raise ValueError(f"DAG edge references unknown rule: '{src}'")
                if dst.lower() not in rule_names:
                    raise ValueError(f"DAG edge references unknown rule: '{dst}'")
        return edges

    @field_validator("rules")
    @classmethod
    def tools_exist_for_rules(cls, rules: list, info) -> list:
        if hasattr(info, 'data') and 'tools' in info.data:
            tool_names = {t.name.lower() for t in info.data["tools"]}
            for rule in rules:
                if rule.tool.lower() not in tool_names:
                    raise ValueError(
                        f"Rule '{rule.name}' references tool '{rule.tool}' "
                        f"which is not in tools list. Available: {tool_names}"
                    )
        return rules

    def get_rule(self, name: str) -> Optional[RuleSpec]:
        return next((r for r in self.rules if r.name == name), None)

    def get_tool(self, name: str) -> Optional[ToolSpec]:
        return next((t for t in self.tools if t.name == name), None)

    def topological_order(self) -> list[str]:
        """Returns rules in DAG execution order (Kahn's algorithm)."""
        from collections import defaultdict, deque
        in_degree = defaultdict(int)
        adjacency = defaultdict(list)
        all_rules  = {r.name for r in self.rules}

        for src, dst in self.dag_edges:
            adjacency[src].append(dst)
            in_degree[dst] += 1

        # Rules with no incoming edges go first
        queue  = deque([r for r in all_rules if in_degree[r] == 0])
        order  = []

        while queue:
            node = queue.popleft()
            order.append(node)
            for neighbor in adjacency[node]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)

        if len(order) != len(all_rules):
            raise ValueError("DAG contains a cycle — cannot determine execution order")

        return order