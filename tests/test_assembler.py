"""
tests/test_assembler.py
=======================
Integration tests for pipeline.assembler — the Jinja2 Snakefile generator.

Tests cover:
  1. render() — pure string output (no filesystem)
  2. assemble() — file writing + dummy input creation
  3. Snakefile content correctness (rule names, wildcards, resources)
  4. config.yaml content correctness
  5. Aggregate rule handling (featurecounts expand())
  6. assemble_from_file() — the JSON → Snakefile convenience path
"""
import re
import pytest
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# 1. render() — pure string output
# ─────────────────────────────────────────────────────────────────────────────

class TestRender:

    @pytest.fixture(scope="class")
    def rendered(self, spec_rna_seq):
        from pipeline.assembler import render
        snakefile, config = render(spec_rna_seq, source_name="test")
        return snakefile, config

    def test_render_returns_two_strings(self, rendered):
        snakefile, config = rendered
        assert isinstance(snakefile, str)
        assert isinstance(config, str)

    def test_snakefile_is_nonempty(self, rendered):
        snakefile, _ = rendered
        assert len(snakefile.strip()) > 0

    def test_config_is_nonempty(self, rendered):
        _, config = rendered
        assert len(config.strip()) > 0

    def test_snakefile_has_rule_all(self, rendered):
        snakefile, _ = rendered
        assert "rule all:" in snakefile

    def test_snakefile_contains_all_rule_names(self, rendered, spec_rna_seq):
        snakefile, _ = rendered
        for rule in spec_rna_seq.rules:
            assert f"rule {rule.name}:" in snakefile, (
                f"Rule '{rule.name}' not found in rendered Snakefile"
            )

    def test_snakefile_has_configfile_directive(self, rendered):
        snakefile, _ = rendered
        assert 'configfile: "config.yaml"' in snakefile

    def test_snakefile_has_samples_line(self, rendered):
        snakefile, _ = rendered
        assert 'samples = config["samples"]' in snakefile

    def test_snakefile_has_threads_for_every_rule(self, rendered, spec_rna_seq):
        snakefile, _ = rendered
        thread_count = snakefile.count("threads:")
        # one per rule (rule all has no threads)
        assert thread_count == len(spec_rna_seq.rules)

    def test_snakefile_has_resources_for_every_rule(self, rendered, spec_rna_seq):
        snakefile, _ = rendered
        resource_count = snakefile.count("resources:")
        assert resource_count == len(spec_rna_seq.rules)

    def test_snakefile_has_log_for_every_rule(self, rendered, spec_rna_seq):
        snakefile, _ = rendered
        log_count = snakefile.count("\n    log:")
        assert log_count == len(spec_rna_seq.rules)

    def test_snakefile_has_shell_or_script_for_every_rule(self, rendered, spec_rna_seq):
        snakefile, _ = rendered
        shell_count  = snakefile.count("\n    shell:")
        script_count = snakefile.count("\n    script:")
        assert (shell_count + script_count) == len(spec_rna_seq.rules)

    def test_source_name_in_snakefile_header(self, rendered):
        snakefile, _ = rendered
        assert "test" in snakefile  # source_name="test" should appear in header comment


# ─────────────────────────────────────────────────────────────────────────────
# 2. Snakefile content correctness
# ─────────────────────────────────────────────────────────────────────────────

class TestSnakefileContent:

    @pytest.fixture(scope="class")
    def snakefile(self, spec_rna_seq):
        from pipeline.assembler import render
        sf, _ = render(spec_rna_seq)
        return sf

    def test_trim_reads_uses_sample_wildcard(self, snakefile):
        """trim_reads is a per-sample rule — must use {sample} wildcard."""
        # Find the trim_reads rule block
        assert re.search(r'rule trim_reads:.*?{sample}', snakefile, re.DOTALL)

    def test_align_star_uses_sample_wildcard(self, snakefile):
        assert re.search(r'rule align_star:.*?{sample}', snakefile, re.DOTALL)

    def test_featurecounts_uses_expand(self, snakefile):
        """featurecounts is aggregate — inputs must be wrapped in expand()."""
        fc_block = _extract_rule_block(snakefile, "featurecounts")
        assert "expand(" in fc_block, (
            "featurecounts rule should use expand() for its inputs"
        )

    def test_featurecounts_no_sample_wildcard_in_output(self, snakefile):
        """Aggregate rule outputs must NOT contain {sample} wildcard."""
        fc_block = _extract_rule_block(snakefile, "featurecounts")
        # Extract output section from the block
        output_match = re.search(r'output:(.*?)(?:params:|threads:|log:)', fc_block, re.DOTALL)
        if output_match:
            output_section = output_match.group(1)
            assert "{sample}" not in output_section

    def test_star_mem_mb_is_large(self, snakefile):
        """STAR needs a lot of RAM — verify mem_mb is written correctly."""
        star_block = _extract_rule_block(snakefile, "align_star")
        assert "40000" in star_block, "STAR rule should have mem_mb=40000"

    def test_log_redirect_in_shell(self, snakefile):
        """Every shell command should redirect stderr to {log}."""
        shell_lines = [l for l in snakefile.splitlines() if '"' in l and "2> {log}" in l]
        assert len(shell_lines) > 0, "No shell commands found with 2> {log} redirect"

    def test_paired_end_inputs_named_r1_r2(self, snakefile):
        """Paired-end FASTQ inputs should be named r1/r2."""
        trim_block = _extract_rule_block(snakefile, "trim_reads")
        assert "r1 =" in trim_block or "r1=" in trim_block
        assert "r2 =" in trim_block or "r2=" in trim_block

    def test_rule_all_targets_featurecounts_output(self, snakefile):
        """rule all must depend on featurecounts outputs (terminal rule)."""
        all_block = _extract_rule_block(snakefile, "all")
        assert "counts/" in all_block


# ─────────────────────────────────────────────────────────────────────────────
# 3. config.yaml content correctness
# ─────────────────────────────────────────────────────────────────────────────

class TestConfigYaml:

    @pytest.fixture(scope="class")
    def config(self, spec_rna_seq):
        from pipeline.assembler import render
        _, cfg = render(spec_rna_seq)
        return cfg

    def test_config_has_samples_section(self, config):
        assert "samples:" in config

    def test_config_lists_all_samples(self, config, spec_rna_seq):
        samples = spec_rna_seq.config_params.get("samples", [])
        for s in samples:
            assert s in config, f"Sample '{s}' missing from config.yaml"

    def test_config_has_genome_reference(self, config):
        """Genome dir should appear in the reference section."""
        assert "star_index" in config or "genome_dir" in config

    def test_config_has_gtf_reference(self, config):
        assert "gtf" in config or "annotation" in config


# ─────────────────────────────────────────────────────────────────────────────
# 4. assemble() — file writing
# ─────────────────────────────────────────────────────────────────────────────

class TestAssemble:

    def test_snakefile_is_written_to_disk(self, assembled_output):
        snakefile, _ = assembled_output
        assert Path(snakefile).exists()

    def test_config_yaml_is_written_to_disk(self, assembled_output):
        _, config = assembled_output
        assert Path(config).exists()

    def test_snakefile_filename(self, assembled_output):
        snakefile, _ = assembled_output
        assert Path(snakefile).name == "Snakefile"

    def test_config_filename(self, assembled_output):
        _, config = assembled_output
        assert Path(config).name == "config.yaml"

    def test_dummy_fastq_files_created(self, assembled_output, spec_rna_seq):
        """assemble() with create_dummy_inputs=True should create .fastq.gz files."""
        snakefile, _ = assembled_output
        out_dir = Path(snakefile).parent
        samples = spec_rna_seq.config_params.get("samples", [])
        for sample in samples:
            r1 = out_dir / "data" / "raw" / f"{sample}_R1.fastq.gz"
            r2 = out_dir / "data" / "raw" / f"{sample}_R2.fastq.gz"
            assert r1.exists(), f"Dummy R1 not created for {sample}"
            assert r2.exists(), f"Dummy R2 not created for {sample}"

    def test_snakefile_content_matches_render(self, assembled_output, spec_rna_seq):
        """Written Snakefile content must match what render() returns."""
        from pipeline.assembler import render
        snakefile_path, _ = assembled_output
        expected, _ = render(spec_rna_seq, source_name="test_fixture")
        actual = Path(snakefile_path).read_text(encoding="utf-8")
        assert actual == expected


# ─────────────────────────────────────────────────────────────────────────────
# 5. assemble_from_file() — JSON → Snakefile path
# ─────────────────────────────────────────────────────────────────────────────

class TestAssembleFromFile:

    def test_loads_and_assembles_fixture(self, tmp_output_dir):
        from pipeline.assembler import assemble_from_file
        fixture = Path(__file__).parent / "fixtures" / "rna_seq_pe.json"
        snakefile, config = assemble_from_file(str(fixture), output_dir=str(tmp_output_dir))
        assert Path(snakefile).exists()
        assert Path(config).exists()

    def test_raises_on_missing_file(self, tmp_output_dir):
        from pipeline.assembler import assemble_from_file
        with pytest.raises(FileNotFoundError):
            assemble_from_file("nonexistent/path.json", output_dir=str(tmp_output_dir))

    def test_raises_on_invalid_json(self, tmp_output_dir, tmp_path):
        from pipeline.assembler import assemble_from_file
        bad_json = tmp_path / "bad.json"
        bad_json.write_text('{"not": "a valid PipelineSpec"}', encoding="utf-8")
        with pytest.raises(Exception):  # ValidationError or similar
            assemble_from_file(str(bad_json), output_dir=str(tmp_output_dir))


# ─────────────────────────────────────────────────────────────────────────────
# 6. Sanitiser unit tests
# ─────────────────────────────────────────────────────────────────────────────

class TestSanitiseShellCmd:
    """
    Tests for _sanitise_shell_cmd — the internal shell command cleaner.
    These are white-box unit tests since the sanitiser is the most fragile part.
    """

    def _sanitise(self, cmd):
        from pipeline.assembler import _sanitise_shell_cmd
        return _sanitise_shell_cmd(cmd)

    def test_strips_trailing_log_redirect(self):
        cmd = "samtools sort -o {output} {input} 2> {log}"
        result = self._sanitise(cmd)
        assert "2> {log}" not in result

    def test_strips_combined_log_redirect(self):
        cmd = "star --runMode genomeGenerate > {log} 2>&1"
        result = self._sanitise(cmd)
        assert "> {log}" not in result

    def test_collapses_wildcard_whitespace(self):
        cmd = "samtools view { sample }.bam"
        result = self._sanitise(cmd)
        assert "{sample}" in result
        assert "{ sample }" not in result

    def test_collapses_double_spaces(self):
        cmd = "fastqc  --outdir  qc/"
        result = self._sanitise(cmd)
        assert "  " not in result

    def test_empty_string_passthrough(self):
        assert self._sanitise("") == ""

    def test_none_passthrough(self):
        # Should handle None gracefully
        from pipeline.assembler import _sanitise_shell_cmd
        result = _sanitise_shell_cmd(None)
        assert result is None or result == ""


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _extract_rule_block(snakefile: str, rule_name: str) -> str:
    """
    Extract the text of a single Snakemake rule from the rendered Snakefile.
    Returns everything from 'rule <name>:' up to the next 'rule ' or end of file.
    """
    pattern = rf"rule {re.escape(rule_name)}:(.*?)(?=\nrule |\Z)"
    match   = re.search(pattern, snakefile, re.DOTALL)
    assert match, f"Rule '{rule_name}' not found in Snakefile"
    return match.group(0)
