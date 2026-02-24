# Workflow Pattern: Bulk RNA-seq Differential Expression Analysis
Category: rna-seq
Analysis type: differential expression, transcriptomics
Keywords: RNA-seq, DESeq2, edgeR, STAR, HISAT2, featureCounts, differential expression, DEG, transcriptome, gene expression

## Purpose
Identify genes that are significantly differentially expressed between two or more conditions
(e.g. treatment vs control, disease vs healthy) from bulk RNA-seq data.

## Pipeline Steps (in order)
1. raw_fastqc         — FastQC QC on raw reads (per sample)
2. trim_reads         — Trimmomatic: adapter removal, quality trimming
3. post_trim_fastqc   — FastQC QC on trimmed reads
4. align_star         — STAR: splice-aware alignment to reference genome
5. sort_index_bam     — samtools sort + index BAM files
6. featurecounts      — featureCounts: gene-level read quantification
7. differential_expr  — DESeq2 (R): DE analysis, generates results table + MA plot + PCA
8. go_enrichment      — clusterProfiler (R): GO and KEGG pathway enrichment
9. multiqc_report     — MultiQC: aggregate all QC logs into HTML report

## Key Decision Points
- STAR vs HISAT2: STAR is faster and more accurate for most use cases; HISAT2 uses less RAM (good for <32GB machines)
- featureCounts vs HTSeq: featureCounts is 10-100x faster; always prefer it
- DESeq2 vs edgeR: DESeq2 recommended for n<20 per group; edgeR equally valid for larger n
- Stranded vs unstranded: check library prep kit; use --stranded 2 for dUTP/reverse-stranded

## File Format Flow
.fastq.gz → [trim] → .fastq.gz → [STAR] → .bam → [sort] → sorted.bam + .bai → [featureCounts] → count_matrix.tsv → [DESeq2] → results.csv

## Required Config Parameters
- samples: list of sample names
- reference: path to reference genome (hg38 / mm10 / etc.)
- annotation: path to GTF file
- metadata: path to sample metadata CSV (must have 'sample' and 'condition' columns)
- adapters: path to adapter FASTA for Trimmomatic

## Resource Requirements (per rule, approximate)
- STAR alignment: 32GB RAM, 8 CPUs, ~45 min per sample
- featureCounts: 4GB RAM, 4 CPUs, ~10 min
- DESeq2: 8GB RAM, 4 CPUs, ~20 min

## Common Errors
- STAR: "genome index not found" — check genomeDir path in config
- featureCounts: low assignment rate (<50%) — check strandedness setting
- DESeq2: "size factors" warning — normalize counts before analysis