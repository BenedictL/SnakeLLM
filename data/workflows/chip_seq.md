# Workflow Pattern: ChIP-seq Histone Modification and Transcription Factor Analysis
Category: chip-seq
Analysis type: histone modification, transcription factor binding, peak calling, epigenomics
Keywords: ChIP-seq, histone, H3K27ac, H3K4me3, transcription factor, TF binding, peak calling, MACS2, ChIP, immunoprecipitation, enrichment, epigenomics

## Purpose
Map genome-wide binding sites of transcription factors or histone modifications
using chromatin immunoprecipitation followed by sequencing (ChIP-seq).
Identifies regulatory elements, promoters, enhancers, and repressed regions.

## Pipeline Steps (in order)
1. raw_fastqc          — FastQC QC on raw reads (input + IP samples)
2. trim_reads          — Trimmomatic: adapter trimming
3. align_bowtie2       — Bowtie2: align to reference genome (single-end or paired-end)
4. filter_bam          — samtools: remove unmapped, low MAPQ (<30), blacklisted reads
5. remove_duplicates   — Picard MarkDuplicates: remove PCR duplicates
6. call_peaks          — MACS2 callpeak: peak calling vs input control
7. filter_blacklist    — bedtools intersect: remove ENCODE blacklist regions
8. create_bigwig       — deeptools bamCompare: IP/input log2 fold-change bigWig tracks
9. peak_annotation     — ChIPseeker (R): annotate peaks to genomic features
10. motif_enrichment   — HOMER findMotifsGenome or MEME-ChIP: TF motif analysis
11. diff_binding       — DiffBind (R): differential binding between conditions
12. multiqc_report     — MultiQC + deeptools plotFingerprint: QC aggregate

## ChIP-seq vs ATAC-seq Key Differences
- ChIP-seq requires INPUT control sample (genomic DNA without antibody pull-down)
- No read shifting required (unlike ATAC-seq)
- MACS2 is used WITH input control: macs2 callpeak -t IP.bam -c input.bam
- Fragment size is ~200bp (mononucleosomal); use MACS2 --nomodel only if sonication QC fails
- Histone marks (broad): H3K27me3, H3K9me3 → use --broad flag in MACS2
- TF binding and active marks (narrow): H3K27ac, H3K4me3, H3K4me1 → default MACS2

## Key Decision Points
- Narrow vs broad peaks: TF/active marks → narrow; repressive marks → --broad
- Paired-end vs single-end: PE recommended; allows fragment size estimation
- Spike-in normalisation: required for quantitative comparison of histone marks between conditions
- IDR: use for TF ChIP replicates; less critical for broad histone marks

## File Format Flow
.fastq.gz → [bowtie2] → .bam → [filter+dedup] → .bam → [MACS2 with input] → .narrowPeak/.broadPeak → [annotate+motif] → results/

## Required Config Parameters
- samples: list of IP sample names
- input_samples: matching input control sample names (dict: ip_sample → input_sample)
- reference: path to reference genome
- blacklist: path to ENCODE blacklist BED file
- effective_genome_size: e.g. 2913022398 for hg38
- peak_type: narrow or broad (depends on antibody target)

## Resource Requirements
- Bowtie2: 8GB RAM, 8 CPUs, ~30 min per sample
- MACS2: 8GB RAM, 4 CPUs, ~15 min
- DiffBind: 16GB RAM, 4 CPUs