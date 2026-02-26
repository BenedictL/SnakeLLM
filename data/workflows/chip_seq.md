# Workflow Pattern: ChIP-seq Chromatin Immunoprecipitation Analysis
Category: chip-seq
Analysis type: protein-DNA binding, epigenomics, transcription factor, histone modification
Keywords: ChIP-seq, peak calling, MACS2, transcription factor, histone modification, IDR, motif, DiffBind, binding

## Purpose
Identify genome-wide binding sites for transcription factors or map histone modifications.
Unlike ATAC-seq, ChIP-seq involves immunoprecipitation of specific proteins cross-linked to DNA or mapping nucleosome bounds.

## Pipeline Steps (in order)
1. raw_fastqc          — FastQC QC on raw reads
2. trim_reads          — Trimmomatic or fastp: adapter trimming
3. align_bowtie2       — Bowtie2: alignment to reference genome 
4. filter_bam          — samtools: remove unmapped, low MAPQ (<20) reads and duplicates
5. remove_duplicates   — Picard MarkDuplicates: remove PCR duplicates
6. call_peaks          — MACS2 callpeak: identify binding peaks (use --broad for histone marks, --narrow for TFs)
7. peak_qc             — deepTools plotFingerprint or plotCoverage: evaluate ChIP-seq quality metrics
8. create_bigwig       — deeptools bamCoverage: generate bigWig tracks for visualization
9. motif_analysis      — HOMER or MEME-ChIP: transcription factor motif enrichment in peaks
10. diff_binding       — DiffBind or DESeq2 on peak counts: differential binding between conditions
11. annotate_peaks     — ChIPseeker (R): annotate peaks to nearest gene, promoter, enhancer
12. multiqc_report     — MultiQC: aggregate all QC

## Critical ChIP-seq Specific Steps (do NOT skip)
- MACS2 requires a control sample (Input DNA or IgG) for accurate peak calling to model background noise.
- Use `--broad` in MACS2 for dispersed marks like H3K36me3, and default (narrow) for sharp TF binding (e.g., CTCF).
- ChIP-seq reads are NOT shifted like ATAC-seq reads. Normal MACS2 model generation is preferred.

## Key Decision Points
- Bowtie2 vs BWA: Bowtie2 preferred for short reads (TF ChIP), BWA-MEM for longer reads.
- IDR filtering: use IDR (Irreproducibility Discovery Rate) if you have biological replicates instead of simply merging peaks.

## File Format Flow
.fastq.gz → [trim] → .fastq.gz → [bowtie2] → .bam → [filter+dedup] → .bam → [MACS2 (with Input)] → .narrowPeak/.broadPeak → [annotate] → peaks_annotated.csv

## Required Config Parameters
- samples: list of sample names (ChIP and Input)
- reference: path to reference genome
- blacklist: path to ENCODE blacklist BED file (removes artifact regions)

## Resource Requirements
- Bowtie2: 8GB RAM, 8 CPUs, ~30 min per sample
- MACS2: 8GB RAM, 4 CPUs, ~15 min
- deeptools bamCoverage: 8GB RAM, 8 CPUs, ~30 min
