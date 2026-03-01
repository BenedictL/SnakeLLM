import os
import shutil
import pytest
import subprocess
from generate_snakefile import generate_pipeline, clean_shell_cmd

def test_clean_shell_cmd():
    """Does the regex cleaner properly fix syntax issues?"""
    # 1. Spaced placeholders
    cmd1 = "cat { input [ 0 ] } > { output }"
    assert clean_shell_cmd(cmd1) == "cat {input[0]} > {output}"
    
    # 2. Bash if/then removal
    cmd2 = "STAR --runThreadN 4 $(if [ -f {input} ]; then echo 'yes'; fi)"
    assert clean_shell_cmd(cmd2) == "STAR --runThreadN 4"
    
    # 3. Log redirects removal
    cmd3 = "bwa mem ref.fa reads.fq > {log} 2>&1"
    assert clean_shell_cmd(cmd3) == "bwa mem ref.fa reads.fq"

def test_snakefile_generate_and_dry_run(tmp_path):
    """Does the generated Snakefile actually pass snakemake --dry-run?"""
    
    valid_spec = {
        "pipeline_type": "unit_test_pipeline",
        "description": "A unit test pipeline",
        "tools": [
            {
                "name": "DummyTool", "version": "1.0", "purpose": "Testing",
                "container": {"registry": "docker.io", "image": "ubuntu", "tag": "latest", "full_uri": "docker.io/ubuntu"}
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
            "samples": ["sample1", "sample2"]
        }
    }
    
    # 1. Generate the pipeline in a temporary pytest directory
    test_out = str(tmp_path / "test_pipeline_out")
    generate_pipeline(valid_spec, test_out, source_name="pytest_mock")
    
    # 2. Assert files were created successfully
    assert os.path.exists(os.path.join(test_out, "Snakefile"))
    assert os.path.exists(os.path.join(test_out, "config.yaml"))
    
    # 3. Create the dummy inputs
    os.makedirs(os.path.join(test_out, "data", "raw"), exist_ok=True)
    with open(os.path.join(test_out, "data", "raw", "input.txt"), "w") as f:
        f.write("dummy")
    with open(os.path.join(test_out, "data", "raw", "intermediate.txt"), "w") as f:
        f.write("dummy")
        
    # 4. Trigger Snakemake dry-run
    # We use subprocess directly to interact with the installed Snakemake binary
    result = subprocess.run(
        ["snakemake", "--dry-run", "--cores", "1"],
        cwd=test_out,
        capture_output=True,
        text=True
    )
    
    # 5. Assert the dry-run completed without Syntax/DAG Errors
    assert result.returncode == 0, f"Snakemake dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    assert "Building DAG of jobs..." in result.stdout
    assert "dummy_rule_two" in result.stdout
