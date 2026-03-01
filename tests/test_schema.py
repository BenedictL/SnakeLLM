import json
import pytest
from pydantic import ValidationError
from core.schema import PipelineSpec, RuleSpec, ToolSpec

def test_pipelinespec_schema_roundtrip():
    """Does a known good JSON payload survive Pydantic parsing and validation?"""
    valid_json = {
        "pipeline_type": "test_pipeline",
        "description": "A unit test pipeline",
        "tools": [
            {
                "name": "DummyTool",
                "version": "1.0",
                "purpose": "Testing",
                "container": {
                    "registry": "docker.io",
                    "image": "ubuntu",
                    "tag": "latest",
                    "full_uri": "docker.io/ubuntu:latest"
                }
            }
        ],
        "rules": [
            {
                "name": "dummy_rule_one",
                "tool": "DummyTool",
                "input": ["data/input.txt"],
                "output": ["data/intermediate.txt"],
                "shell_cmd": "cat {input} > {output}"
            },
            {
                "name": "dummy_rule_two",
                "tool": "DummyTool",
                "input": ["data/intermediate.txt"],
                "output": ["data/output.txt"],
                "shell_cmd": "cat {input} > {output}"
            }
        ],
        "dag_edges": [
            ["dummy_rule_one", "dummy_rule_two"]
        ],
        "config_params": {
            "samples": ["sample1"]
        }
    }
    
    # Validation should succeed
    spec = PipelineSpec(**valid_json)
    assert len(spec.rules) == 2
    assert spec.rules[0].name == "dummy_rule_one"
    
def test_schema_cross_validation_catches_missing_tool():
    """Does the @model_validator catch rules mapped to non-existent tools?"""
    invalid_json = {
        "pipeline_type": "test",
        "description": "test",
        "tools": [
            {
                "name": "RealTool",
                "version": "1.0",
                "purpose": "Testing",
                "container": {"registry": "docker.io", "image": "ubuntu", "tag": "latest", "full_uri": "docker.io/ubuntu"}
            }
        ],
        "rules": [
            {
                "name": "bad_rule",
                "tool": "FakeTool",  # This tool doesn't exist in the list above
                "input": [],
                "output": ["out.txt"],
                "shell_cmd": "echo"
            }
        ],
        "dag_edges": [],
        "config_params": {}
    }
    
    with pytest.raises(ValidationError) as excinfo:
        PipelineSpec(**invalid_json)
    
    # Asserting the specific cross-field custom error is triggered
    assert "references tool 'FakeTool' which is not in tools list" in str(excinfo.value)

def test_dag_cycle_detection():
    """Does topological_order() properly raise an error on cyclic DAGs?"""
    cyclic_json = {
        "pipeline_type": "cycle_test",
        "description": "Tests cycle breaking",
        "tools": [
            {
                "name": "T", "version": "1", "purpose": "T",
                "container": {"registry": "r", "image": "i", "tag": "t", "full_uri": "r/i:t"}
            }
        ],
        "rules": [
            {"name": "rule_a", "tool": "T", "input": [], "output": ["a.txt"], "shell_cmd": "echo"},
            {"name": "rule_b", "tool": "T", "input": [], "output": ["b.txt"], "shell_cmd": "echo"},
            {"name": "rule_c", "tool": "T", "input": [], "output": ["c.txt"], "shell_cmd": "echo"}
        ],
        "dag_edges": [
            ["rule_a", "rule_b"],
            ["rule_b", "rule_c"],
            ["rule_c", "rule_a"]  # This makes it a circle A -> B -> C -> A
        ],
        "config_params": {}
    }
    
    spec = PipelineSpec(**cyclic_json)
    
    with pytest.raises(ValueError) as excinfo:
        spec.topological_order()
        
    assert "DAG contains a cycle — cannot determine execution order" in str(excinfo.value)
