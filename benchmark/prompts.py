"""
benchmark/prompts.py
====================
18 benchmark prompts: 6 pipeline types × 3 difficulty levels.

Each prompt has:
  - prompt_id   : unique identifier
  - pipeline_type: which workflow doc it maps to
  - difficulty  : simple | medium | hard
  - prompt      : natural language input to the LLM
  - expected    : ground-truth values for metric evaluation

Difficulty levels:
  simple  — single tool chain, standard terminology, no ambiguity
  medium  — multi-step, some tool name ambiguity, paired-end or multi-sample
  hard    — branching DAG, aggregate steps, domain-specific terminology,
             tool name disambiguation required
"""

from dataclasses import dataclass, field
from typing import Literal

Difficulty = Literal["simple", "medium", "hard"]


@dataclass
class BenchmarkPrompt:
    prompt_id:     str
    pipeline_type: str
    difficulty:    Difficulty
    prompt:        str
    expected_tools:      list[str]   # canonical tool names (lowercase)
    expected_rule_count: int         # minimum expected rules
    expected_dag_edges:  list[tuple[str, str]]  # (src, dst) — subset check
    expected_containers: list[str]   # substrings that must appear in container URIs


PROMPTS: list[BenchmarkPrompt] = [

    # ─────────────────────────────────────────────────────────────────────────
    # RNA-SEQ DE
    # ─────────────────────────────────────────────────────────────────────────
    BenchmarkPrompt(
        prompt_id="rna_seq_simple",
        pipeline_type="rna_seq_de",
        difficulty="simple",
        prompt="Align RNA-seq reads to a reference genome using STAR and quantify gene expression with featureCounts.",
        expected_tools=["star", "featurecounts"],
        expected_rule_count=2,
        expected_dag_edges=[("align_star", "featurecounts")],
        expected_containers=["star", "subread"],
    ),
    BenchmarkPrompt(
        prompt_id="rna_seq_medium",
        pipeline_type="rna_seq_de",
        difficulty="medium",
        prompt="Run a paired-end RNA-seq pipeline: quality control with FastQC, adapter trimming with Trimmomatic, alignment with STAR, and gene quantification with featureCounts.",
        expected_tools=["fastqc", "trimmomatic", "star", "featurecounts"],
        expected_rule_count=4,
        expected_dag_edges=[
            ("trim_reads", "align_star"),
            ("align_star", "featurecounts"),
        ],
        expected_containers=["fastqc", "trimmomatic", "star", "subread"],
    ),
    BenchmarkPrompt(
        prompt_id="rna_seq_hard",
        pipeline_type="rna_seq_de",
        difficulty="hard",
        prompt=(
            "Build a complete differential expression pipeline for paired-end RNA-seq: "
            "FastQC on raw reads, Trimmomatic adapter trimming, STAR genome alignment, "
            "samtools BAM sorting and indexing, featureCounts quantification, and DESeq2 "
            "differential expression analysis. Use Wald test with BH correction, padj < 0.05."
        ),
        expected_tools=["fastqc", "trimmomatic", "star", "samtools", "featurecounts", "deseq2"],
        expected_rule_count=6,
        expected_dag_edges=[
            ("trim_reads", "align_star"),
            ("align_star", "sort_index_bam"),
            ("sort_index_bam", "featurecounts"),
            ("featurecounts", "deseq2"),
        ],
        expected_containers=["fastqc", "trimmomatic", "star", "samtools", "subread", "deseq2"],
    ),

    # ─────────────────────────────────────────────────────────────────────────
    # ATAC-SEQ
    # ─────────────────────────────────────────────────────────────────────────
    BenchmarkPrompt(
        prompt_id="atac_seq_simple",
        pipeline_type="atac_seq",
        difficulty="simple",
        prompt="Call peaks from ATAC-seq BAM files using MACS2.",
        expected_tools=["macs2"],
        expected_rule_count=1,
        expected_dag_edges=[],
        expected_containers=["macs2"],
    ),
    BenchmarkPrompt(
        prompt_id="atac_seq_medium",
        pipeline_type="atac_seq",
        difficulty="medium",
        prompt=(
            "ATAC-seq pipeline: trim adapters with Trimmomatic, align to genome with Bowtie2, "
            "remove duplicates with Picard MarkDuplicates, and call peaks with MACS2."
        ),
        expected_tools=["trimmomatic", "bowtie2", "picard", "macs2"],
        expected_rule_count=4,
        expected_dag_edges=[
            ("trim_reads", "align_bowtie2"),
            ("align_bowtie2", "mark_duplicates"),
            ("mark_duplicates", "call_peaks"),
        ],
        expected_containers=["trimmomatic", "bowtie2", "picard", "macs2"],
    ),
    BenchmarkPrompt(
        prompt_id="atac_seq_hard",
        pipeline_type="atac_seq",
        difficulty="hard",
        prompt=(
            "Full ATAC-seq pipeline with IDR: FastQC, Trimmomatic trimming, Bowtie2 alignment, "
            "samtools filtering (remove chrM, unmapped, low MAPQ), Picard MarkDuplicates, "
            "MACS2 peak calling per replicate, IDR analysis across replicates, "
            "and MultiQC aggregated QC report."
        ),
        expected_tools=["fastqc", "trimmomatic", "bowtie2", "samtools", "picard", "macs2", "multiqc"],
        expected_rule_count=7,
        expected_dag_edges=[
            ("trim_reads", "align_bowtie2"),
            ("align_bowtie2", "filter_bam"),
            ("filter_bam", "mark_duplicates"),
            ("mark_duplicates", "call_peaks"),
        ],
        expected_containers=["fastqc", "trimmomatic", "bowtie2", "samtools", "picard", "macs2", "multiqc"],
    ),

    # ─────────────────────────────────────────────────────────────────────────
    # WGS VARIANT CALLING
    # ─────────────────────────────────────────────────────────────────────────
    BenchmarkPrompt(
        prompt_id="wgs_simple",
        pipeline_type="wgs_variant_calling",
        difficulty="simple",
        prompt="Call germline variants from WGS data using GATK HaplotypeCaller.",
        expected_tools=["gatk"],
        expected_rule_count=1,
        expected_dag_edges=[],
        expected_containers=["gatk"],
    ),
    BenchmarkPrompt(
        prompt_id="wgs_medium",
        pipeline_type="wgs_variant_calling",
        difficulty="medium",
        prompt=(
            "WGS variant calling pipeline: align reads with BWA-MEM2, sort and index with samtools, "
            "mark duplicates with Picard, and call variants with GATK HaplotypeCaller in GVCF mode."
        ),
        expected_tools=["bwa-mem2", "samtools", "picard", "gatk"],
        expected_rule_count=4,
        expected_dag_edges=[
            ("align_bwa", "sort_index_bam"),
            ("sort_index_bam", "mark_duplicates"),
            ("mark_duplicates", "haplotype_caller"),
        ],
        expected_containers=["bwamem2", "samtools", "picard", "gatk"],
    ),
    BenchmarkPrompt(
        prompt_id="wgs_hard",
        pipeline_type="wgs_variant_calling",
        difficulty="hard",
        prompt=(
            "GATK best practices WGS pipeline: FastQC, BWA-MEM2 alignment, samtools sort/index, "
            "Picard MarkDuplicates, GATK BaseRecalibrator + ApplyBQSR, HaplotypeCaller per-sample GVCF, "
            "GenomicsDBImport joint genotyping, GenotypeGVCFs, VQSR filtering, and bcftools stats."
        ),
        expected_tools=["fastqc", "bwa-mem2", "samtools", "picard", "gatk", "bcftools"],
        expected_rule_count=8,
        expected_dag_edges=[
            ("align_bwa", "sort_index_bam"),
            ("sort_index_bam", "mark_duplicates"),
            ("mark_duplicates", "base_recalibration"),
            ("base_recalibration", "haplotype_caller"),
            ("haplotype_caller", "joint_genotyping"),
        ],
        expected_containers=["bwamem2", "samtools", "picard", "gatk", "bcftools"],
    ),

    # ─────────────────────────────────────────────────────────────────────────
    # SCRNA-SEQ
    # ─────────────────────────────────────────────────────────────────────────
    BenchmarkPrompt(
        prompt_id="scrna_simple",
        pipeline_type="scrna_seq",
        difficulty="simple",
        prompt="Quantify single-cell RNA-seq data from 10x Genomics using Cell Ranger count.",
        expected_tools=["cellranger"],
        expected_rule_count=1,
        expected_dag_edges=[],
        expected_containers=["cellranger"],
    ),
    BenchmarkPrompt(
        prompt_id="scrna_medium",
        pipeline_type="scrna_seq",
        difficulty="medium",
        prompt=(
            "scRNA-seq preprocessing pipeline: run Cell Ranger count for each sample, "
            "then aggregate with Cell Ranger aggr for multi-sample analysis."
        ),
        expected_tools=["cellranger"],
        expected_rule_count=2,
        expected_dag_edges=[("cellranger_count", "cellranger_aggr")],
        expected_containers=["cellranger"],
    ),
    BenchmarkPrompt(
        prompt_id="scrna_hard",
        pipeline_type="scrna_seq",
        difficulty="hard",
        prompt=(
            "Full scRNA-seq pipeline: FastQC on raw reads, Cell Ranger count per sample, "
            "Cell Ranger aggr for multi-sample integration, then Seurat clustering analysis "
            "with dimensionality reduction (PCA, UMAP), marker gene identification, "
            "and cell type annotation. Filter cells with >200 genes and <20% mitochondrial reads."
        ),
        expected_tools=["fastqc", "cellranger", "seurat"],
        expected_rule_count=4,
        expected_dag_edges=[
            ("cellranger_count", "cellranger_aggr"),
            ("cellranger_aggr", "seurat_analysis"),
        ],
        expected_containers=["fastqc", "cellranger", "seurat"],
    ),

    # ─────────────────────────────────────────────────────────────────────────
    # CHIP-SEQ
    # ─────────────────────────────────────────────────────────────────────────
    BenchmarkPrompt(
        prompt_id="chip_seq_simple",
        pipeline_type="chip_seq",
        difficulty="simple",
        prompt="Call ChIP-seq peaks using MACS2 with an input control.",
        expected_tools=["macs2"],
        expected_rule_count=1,
        expected_dag_edges=[],
        expected_containers=["macs2"],
    ),
    BenchmarkPrompt(
        prompt_id="chip_seq_medium",
        pipeline_type="chip_seq",
        difficulty="medium",
        prompt=(
            "ChIP-seq pipeline: trim adapters with Trimmomatic, align with Bowtie2, "
            "remove duplicates with Picard, and call peaks with MACS2 using matched input control."
        ),
        expected_tools=["trimmomatic", "bowtie2", "picard", "macs2"],
        expected_rule_count=4,
        expected_dag_edges=[
            ("trim_reads", "align_bowtie2"),
            ("align_bowtie2", "mark_duplicates"),
            ("mark_duplicates", "call_peaks"),
        ],
        expected_containers=["trimmomatic", "bowtie2", "picard", "macs2"],
    ),
    BenchmarkPrompt(
        prompt_id="chip_seq_hard",
        pipeline_type="chip_seq",
        difficulty="hard",
        prompt=(
            "Complete H3K27ac ChIP-seq pipeline: FastQC, Trimmomatic, Bowtie2 alignment, "
            "samtools filtering and blacklist removal, Picard MarkDuplicates, "
            "MACS2 peak calling with input control, deepTools bamCoverage for bigWig tracks, "
            "deepTools computeMatrix and plotHeatmap for visualisation, and MultiQC report."
        ),
        expected_tools=["fastqc", "trimmomatic", "bowtie2", "samtools", "picard", "macs2", "deeptools", "multiqc"],
        expected_rule_count=8,
        expected_dag_edges=[
            ("trim_reads", "align_bowtie2"),
            ("align_bowtie2", "filter_bam"),
            ("filter_bam", "mark_duplicates"),
            ("mark_duplicates", "call_peaks"),
            ("mark_duplicates", "bam_coverage"),
        ],
        expected_containers=["trimmomatic", "bowtie2", "samtools", "picard", "macs2", "deeptools", "multiqc"],
    ),

    # ─────────────────────────────────────────────────────────────────────────
    # METHYLATION / WGBS
    # ─────────────────────────────────────────────────────────────────────────
    BenchmarkPrompt(
        prompt_id="methylation_simple",
        pipeline_type="methylation_wgbs",
        difficulty="simple",
        prompt="Align WGBS reads and extract methylation calls using Bismark.",
        expected_tools=["bismark"],
        expected_rule_count=2,
        expected_dag_edges=[("bismark_align", "bismark_extract")],
        expected_containers=["bismark"],
    ),
    BenchmarkPrompt(
        prompt_id="methylation_medium",
        pipeline_type="methylation_wgbs",
        difficulty="medium",
        prompt=(
            "WGBS methylation pipeline: trim adapters with Trim Galore, "
            "align bisulfite-converted reads with Bismark, deduplicate, "
            "and extract CpG methylation calls."
        ),
        expected_tools=["trim-galore", "bismark"],
        expected_rule_count=4,
        expected_dag_edges=[
            ("trim_reads", "bismark_align"),
            ("bismark_align", "bismark_deduplicate"),
            ("bismark_deduplicate", "bismark_extract"),
        ],
        expected_containers=["trim-galore", "bismark"],
    ),
    BenchmarkPrompt(
        prompt_id="methylation_hard",
        pipeline_type="methylation_wgbs",
        difficulty="hard",
        prompt=(
            "Complete WGBS pipeline: FastQC, Trim Galore adapter and quality trimming, "
            "Bismark genome preparation, bisulfite alignment, deduplication, "
            "methylation extraction for CpG/CHG/CHH contexts, "
            "MethylKit differential methylation analysis across conditions, "
            "and MultiQC report aggregating FastQC, Trimming, and Bismark metrics."
        ),
        expected_tools=["fastqc", "trim-galore", "bismark", "methylkit", "multiqc"],
        expected_rule_count=7,
        expected_dag_edges=[
            ("trim_reads", "bismark_align"),
            ("bismark_align", "bismark_deduplicate"),
            ("bismark_deduplicate", "bismark_extract"),
            ("bismark_extract", "differential_methylation"),
        ],
        expected_containers=["fastqc", "trim-galore", "bismark", "methylkit", "multiqc"],
    ),
]

# Lookup helpers
PROMPT_BY_ID = {p.prompt_id: p for p in PROMPTS}

def get_prompts_by_difficulty(difficulty: Difficulty) -> list[BenchmarkPrompt]:
    return [p for p in PROMPTS if p.difficulty == difficulty]

def get_prompts_by_pipeline(pipeline_type: str) -> list[BenchmarkPrompt]:
    return [p for p in PROMPTS if p.pipeline_type == pipeline_type]
