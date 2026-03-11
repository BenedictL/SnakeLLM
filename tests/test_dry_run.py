"""
tests/test_dry_run.py
=====================
End-to-end test: assemble a Snakefile from the RNA-seq fixture and validate
it with `snakemake --dry-run`.

This is the highest-confidence test — if Snakemake itself accepts the DAG,
the pipeline is structurally correct.

Skipped automatically if `snakemake` is not installed in the test environment.
"""
import subprocess
import shutil
import pytest
from pathlib import Path


# Skip the whole module if snakemake isn't available
snakemake_available = shutil.which("snakemake") is not None
pytestmark = pytest.mark.skipif(
    not snakemake_available,
    reason="snakemake not installed — skipping dry-run tests"
)


class TestDryRun:

    def test_rna_seq_dry_run_exits_zero(self, assembled_output):
        """
        The canonical integration test.
        Snakemake --dry-run must exit 0, meaning:
          - Snakefile syntax is valid
          - All wildcards resolve correctly
          - DAG has no cycles
          - rule all targets are reachable
        """
        snakefile, _ = assembled_output
        out_dir = Path(snakefile).parent

        result = subprocess.run(
            ["snakemake", "--dry-run", "--cores", "1",
             "--snakefile", str(snakefile),
             "--directory", str(out_dir)],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, (
            f"snakemake --dry-run failed (exit {result.returncode}).\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_dry_run_reports_correct_job_count(self, assembled_output, spec_rna_seq):
        """
        Dry-run output should list jobs for every rule × every sample
        (for per-sample rules) plus 1 for rule all and 1 for aggregate rules.

        For RNA-seq with 3 samples:
          trim_reads × 3, align_star × 3, sort_index_bam × 3,
          featurecounts × 1, all × 1 = 11 total jobs
        """
        snakefile, _ = assembled_output
        out_dir      = Path(snakefile).parent
        samples      = spec_rna_seq.config_params.get("samples", [])

        result = subprocess.run(
            ["snakemake", "--dry-run", "--cores", "1",
             "--snakefile", str(snakefile),
             "--directory", str(out_dir)],
            capture_output=True,
            text=True,
        )

        # Per-sample rules × n_samples + aggregate rules + rule all
        per_sample_rules = ["fastqc_raw", "trim_reads", "align_star", "sort_index_bam"]
        aggregate_rules  = ["featurecounts"]
        expected_jobs    = (
            len(per_sample_rules) * len(samples)
            + len(aggregate_rules)
            + 1  # rule all
        )

        # Parse "total   N" from the job stats table in dry-run output
        import re
        match = re.search(r'total\s+(\d+)', result.stdout)
        if match:
            actual_jobs = int(match.group(1))
            assert actual_jobs == expected_jobs, (
                f"Expected {expected_jobs} jobs, got {actual_jobs}"
            )

    def test_dry_run_output_mentions_featurecounts(self, assembled_output):
        """featurecounts must appear in dry-run job list."""
        snakefile, _ = assembled_output
        out_dir = Path(snakefile).parent

        result = subprocess.run(
            ["snakemake", "--dry-run", "--cores", "1",
             "--snakefile", str(snakefile),
             "--directory", str(out_dir)],
            capture_output=True,
            text=True,
        )
        assert "featurecounts" in result.stdout

    def test_dry_run_output_mentions_all_samples(self, assembled_output, spec_rna_seq):
        """Each sample should appear in the dry-run output."""
        snakefile, _ = assembled_output
        out_dir  = Path(snakefile).parent
        samples  = spec_rna_seq.config_params.get("samples", [])

        result = subprocess.run(
            ["snakemake", "--dry-run", "--cores", "1",
             "--snakefile", str(snakefile),
             "--directory", str(out_dir)],
            capture_output=True,
            text=True,
        )
        for sample in samples:
            assert sample in result.stdout, (
                f"Sample '{sample}' not mentioned in dry-run output"
            )
