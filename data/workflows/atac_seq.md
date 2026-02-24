# Workflow Pattern: ATAC-seq Chromatin Accessibility Analysis
Category: atac-seq
Analysis type: chromatin accessibility, open chromatin, epigenomics
Keywords: ATAC-seq, peak calling, MACS2, chromatin, open chromatin, accessibility, transposase, nucleosome, IDR, motif

## Purpose
Map genome-wide open chromatin regions (accessible DNA) using the ATAC-seq assay.
Identifies regulatory elements, enhancers, and promoters that are accessible for transcription factor binding.

## Pipeline Steps (in order)
1. raw_fastqc          — FastQC QC on raw reads
2. trim_reads          — Trimmomatic or fastp: adapter trimming (use Nextera adapters for ATAC-seq)
3. align_bowtie2       — Bowtie2: alignment to reference genome (NOT STAR — ATAC-seq is not RNA)
4. filter_bam          — samtools: remove unmapped, low MAPQ (<30), mitochondrial reads
5. remove_duplicates   — Picard MarkDuplicates: remove PCR duplicates
6. shift_reads         — deeptools alignmentSieve: shift reads +4/-5 bp (Tn5 transposase correction)
7. call_peaks          — MACS2 callpeak: identify open chromatin peaks (--nomodel --shift -100 --extsize 200)
8. peak_qc             — ataqv: ATAC-seq specific QC (TSS enrichment, fragment size distribution)
9. create_bigwig       — deeptools bamCoverage: generate bigWig tracks for visualization
10. motif_analysis     — HOMER or MEME-ChIP: transcription factor motif enrichment in peaks
11. diff_accessibility — DiffBind or DESeq2 on peak counts: differential accessibility between conditions
12. annotate_peaks     — ChIPseeker (R): annotate peaks to nearest gene, promoter, enhancer
13. multiqc_report     — MultiQC: aggregate all QC

## Critical ATAC-seq Specific Steps (do NOT skip)
- Read shifting (+4 bp on + strand, -5 bp on - strand) is mandatory — corrects for Tn5 insertion bias
- Remove mitochondrial reads — MT chromosome is highly accessible and inflates background
- Fragment size distribution should show nucleosomal banding (~200bp, ~400bp, ~600bp peaks)
- TSS enrichment score should be >7 for high-quality ATAC-seq

## Key Decision Points
- Bowtie2 vs BWA: Bowtie2 preferred for ATAC-seq (short reads, no splicing)
- MACS2 vs MACS3: MACS3 is the updated version; both valid; use --nomodel for ATAC-seq
- IDR filtering: use IDR (Irreproducibility Discovery Rate) if you have biological replicates

## File Format Flow
.fastq.gz → [trim] → .fastq.gz → [bowtie2] → .bam → [filter+dedup+shift] → .bam → [MACS2] → .narrowPeak → [annotate] → peaks_annotated.csv

## Required Config Parameters
- samples: list of sample names
- reference: path to reference genome
- blacklist: path to ENCODE blacklist BED file (removes artifact regions)
- tss_bed: path to TSS annotation BED file (for TSS enrichment QC)
- effective_genome_size: e.g. 2913022398 for hg38 (needed by MACS2 and deeptools)

## Resource Requirements
- Bowtie2: 8GB RAM, 8 CPUs, ~30 min per sample
- MACS2: 16GB RAM, 4 CPUs, ~20 min
- deeptools bamCoverage: 8GB RAM, 8 CPUs, ~30 min