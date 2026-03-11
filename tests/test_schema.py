"""
tests/test_schema.py
====================
Unit tests for core.schema — PipelineSpec, RuleSpec, ResourceSpec validation.

These tests do NOT call the LLM or touch the filesystem (beyond loading the
fixture JSON). They verify that:
  1. Valid JSON round-trips cleanly through Pydantic
  2. Required fields are present and correctly typed
  3. Derived methods (topological_order, dag_edges) behave correctly
  4. Invalid specs are rejected with clear errors
"""
import json
import pytest
from pathlib import Path
from pydantic import ValidationError

FIXTURE_DIR = Path(__file__).parent / "fixtures"


# ─────────────────────────────────────────────────────────────────────────────
# 1. Basic loading & field presence
# ─────────────────────────────────────────────────────────────────────────────

class TestPipelineSpecLoading:

    def test_loads_from_fixture_json(self, spec_rna_seq):
        """PipelineSpec loads without errors from the RNA-seq fixture."""
        assert spec_rna_seq is not None

    def test_pipeline_type_is_string(self, spec_rna_seq):
        assert isinstance(spec_rna_seq.pipeline_type, str)
        assert len(spec_rna_seq.pipeline_type) > 0

    def test_description_is_string(self, spec_rna_seq):
        assert isinstance(spec_rna_seq.description, str)

    def test_rules_is_nonempty_list(self, spec_rna_seq):
        assert isinstance(spec_rna_seq.rules, list)
        assert len(spec_rna_seq.rules) > 0

    def test_tools_is_nonempty_list(self, spec_rna_seq):
        assert isinstance(spec_rna_seq.tools, list)
        assert len(spec_rna_seq.tools) > 0

    def test_config_params_is_dict(self, spec_rna_seq):
        assert isinstance(spec_rna_seq.config_params, dict)

    def test_samples_in_config_params(self, spec_rna_seq):
        assert "samples" in spec_rna_seq.config_params
        samples = spec_rna_seq.config_params["samples"]
        assert isinstance(samples, list)
        assert len(samples) > 0

    def test_json_round_trip(self, spec_rna_seq):
        """Serialise → deserialise must produce an equal spec."""
        from core.schema import PipelineSpec
        serialised   = spec_rna_seq.model_dump_json()
        deserialised = PipelineSpec.model_validate_json(serialised)
        assert deserialised.pipeline_type == spec_rna_seq.pipeline_type
        assert len(deserialised.rules)    == len(spec_rna_seq.rules)


# ─────────────────────────────────────────────────────────────────────────────
# 2. RuleSpec field validation
# ─────────────────────────────────────────────────────────────────────────────

class TestRuleSpec:

    def test_all_rules_have_names(self, spec_rna_seq):
        for rule in spec_rna_seq.rules:
            assert isinstance(rule.name, str)
            assert len(rule.name) > 0

    def test_rule_names_are_unique(self, spec_rna_seq):
        names = [r.name for r in spec_rna_seq.rules]
        assert len(names) == len(set(names)), "Duplicate rule names found"

    def test_all_rules_have_inputs(self, spec_rna_seq):
        for rule in spec_rna_seq.rules:
            assert isinstance(rule.input, list)
            assert len(rule.input) > 0, f"Rule {rule.name} has empty input"

    def test_all_rules_have_outputs(self, spec_rna_seq):
        for rule in spec_rna_seq.rules:
            assert isinstance(rule.output, list)
            assert len(rule.output) > 0, f"Rule {rule.name} has empty output"

    def test_threads_is_positive_int(self, spec_rna_seq):
        for rule in spec_rna_seq.rules:
            if rule.threads is not None:
                assert isinstance(rule.threads, int)
                assert rule.threads > 0, f"Rule {rule.name} has non-positive threads"

    def test_resources_mem_mb_positive(self, spec_rna_seq):
        for rule in spec_rna_seq.rules:
            assert rule.resources.mem_mb > 0, f"Rule {rule.name} has zero mem_mb"

    def test_specific_rules_present(self, spec_rna_seq):
        """RNA-seq pipeline must contain these exact rules."""
        names = {r.name for r in spec_rna_seq.rules}
        expected = {"trim_reads", "align_star", "sort_index_bam", "featurecounts"}
        assert expected.issubset(names), f"Missing rules: {expected - names}"


# ─────────────────────────────────────────────────────────────────────────────
# 3. DAG structure & topological order
# ─────────────────────────────────────────────────────────────────────────────

class TestDAGStructure:

    def test_dag_edges_is_list(self, spec_rna_seq):
        assert isinstance(spec_rna_seq.dag_edges, list)

    def test_dag_edges_are_pairs(self, spec_rna_seq):
        for edge in spec_rna_seq.dag_edges:
            assert len(edge) == 2, f"Edge {edge} is not a (src, dst) pair"

    def test_dag_edge_nodes_exist_as_rules(self, spec_rna_seq):
        """Every node mentioned in dag_edges must correspond to a rule name."""
        rule_names = {r.name for r in spec_rna_seq.rules}
        for src, dst in spec_rna_seq.dag_edges:
            assert src in rule_names, f"DAG source '{src}' has no matching rule"
            assert dst in rule_names, f"DAG dest '{dst}' has no matching rule"

    def test_topological_order_returns_all_rules(self, spec_rna_seq):
        order = spec_rna_seq.topological_order()
        rule_names = {r.name for r in spec_rna_seq.rules}
        assert set(order) == rule_names

    def test_topological_order_respects_edges(self, spec_rna_seq):
        """For every edge (src→dst), src must appear before dst in the order."""
        order = spec_rna_seq.topological_order()
        position = {name: i for i, name in enumerate(order)}
        for src, dst in spec_rna_seq.dag_edges:
            assert position[src] < position[dst], (
                f"Topological order violated: {src} (pos {position[src]}) "
                f"should come before {dst} (pos {position[dst]})"
            )

    def test_trim_reads_before_align_star(self, spec_rna_seq):
        order = spec_rna_seq.topological_order()
        assert order.index("trim_reads") < order.index("align_star")

    def test_align_star_before_sort_index_bam(self, spec_rna_seq):
        order = spec_rna_seq.topological_order()
        assert order.index("align_star") < order.index("sort_index_bam")

    def test_sort_index_bam_before_featurecounts(self, spec_rna_seq):
        order = spec_rna_seq.topological_order()
        assert order.index("sort_index_bam") < order.index("featurecounts")

    def test_featurecounts_is_terminal(self, spec_rna_seq):
        """featurecounts is a sink node — nothing should depend on it."""
        sources = {src for src, _ in spec_rna_seq.dag_edges}
        assert "featurecounts" not in sources


# ─────────────────────────────────────────────────────────────────────────────
# 4. Rejection of invalid specs
# ─────────────────────────────────────────────────────────────────────────────

class TestSchemaValidationErrors:

    def test_rejects_rules_with_min_length(self, rna_seq_json):
        """PipelineSpec.rules has min_length=1 — empty list must raise."""
        from core.schema import PipelineSpec
        bad = {**rna_seq_json, "rules": []}
        with pytest.raises(ValidationError):
            PipelineSpec.model_validate(bad)

    def test_rejects_tools_with_min_length(self, rna_seq_json):
        """PipelineSpec.tools has min_length=1 — empty list must raise."""
        from core.schema import PipelineSpec
        bad = {**rna_seq_json, "tools": []}
        with pytest.raises(ValidationError):
            PipelineSpec.model_validate(bad)

    def test_rejects_invalid_rule_name(self, rna_seq_json):
        """RuleSpec.name must be snake_case — CamelCase should raise."""
        from core.schema import PipelineSpec
        import copy
        bad = copy.deepcopy(rna_seq_json)
        bad["rules"][0]["name"] = "TrimReads"  # CamelCase — fails snake_case validator
        with pytest.raises(ValidationError):
            PipelineSpec.model_validate(bad)

    def test_rejects_rule_referencing_unknown_tool(self, rna_seq_json):
        """model_validator checks rule.tool is in tools list."""
        from core.schema import PipelineSpec
        import copy
        bad = copy.deepcopy(rna_seq_json)
        bad["rules"][0]["tool"] = "nonexistent_tool"
        with pytest.raises((ValidationError, ValueError)):
            PipelineSpec.model_validate(bad)

    def test_rejects_dag_edge_referencing_unknown_rule(self, rna_seq_json):
        """model_validator checks all dag_edges reference existing rules."""
        from core.schema import PipelineSpec
        import copy
        bad = copy.deepcopy(rna_seq_json)
        bad["dag_edges"].append(["trim_reads", "nonexistent_rule"])
        with pytest.raises((ValidationError, ValueError)):
            PipelineSpec.model_validate(bad)

    def test_rejects_mem_mb_below_minimum(self, rna_seq_json):
        """ResourceSpec.mem_mb has ge=256 — values below must raise."""
        from core.schema import PipelineSpec
        import copy
        bad = copy.deepcopy(rna_seq_json)
        bad["rules"][0]["resources"]["mem_mb"] = 100  # below ge=256
        with pytest.raises(ValidationError):
            PipelineSpec.model_validate(bad)

    def test_rejects_time_min_below_minimum(self, rna_seq_json):
        """ResourceSpec.time_min has ge=1 — zero must raise."""
        from core.schema import PipelineSpec
        import copy
        bad = copy.deepcopy(rna_seq_json)
        bad["rules"][0]["resources"]["time_min"] = 0  # below ge=1
        with pytest.raises(ValidationError):
            PipelineSpec.model_validate(bad)
