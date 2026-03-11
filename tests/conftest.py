"""
tests/conftest.py
=================
Shared pytest fixtures for SnakeLLM integration tests.

Provides:
  - spec_rna_seq     : validated PipelineSpec loaded from fixture JSON
  - tmp_output_dir   : fresh temp directory per test
  - assembled_output : (snakefile_path, config_path) written to tmp dir
"""
import json
import pytest
from pathlib import Path

# ── Fixture directory ─────────────────────────────────────────────────────────
FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def rna_seq_json() -> dict:
    """Raw dict from the RNA-seq fixture JSON — not yet validated."""
    path = FIXTURE_DIR / "rna_seq_pe.json"
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="session")
def spec_rna_seq(rna_seq_json):
    """
    Validated PipelineSpec for the paired-end RNA-seq pipeline.

    Scope=session: Pydantic validation runs once for the whole test session.
    This fixture is the shared ground-truth object for all schema/assembler tests.
    """
    from core.schema import PipelineSpec
    return PipelineSpec.model_validate(rna_seq_json)


@pytest.fixture()
def tmp_output_dir(tmp_path) -> Path:
    """
    Fresh temporary directory for each test.
    pytest's built-in tmp_path gives an isolated dir per test invocation.
    """
    return tmp_path / "output"


@pytest.fixture()
def assembled_output(spec_rna_seq, tmp_output_dir):
    """
    Runs the assembler once and returns (snakefile_path, config_path).
    Use this in any test that needs the rendered files on disk.
    """
    from pipeline.assembler import assemble
    snakefile, config = assemble(
        spec_rna_seq,
        output_dir=str(tmp_output_dir),
        source_name="test_fixture",
        create_dummy_inputs=True,
    )
    return snakefile, config
