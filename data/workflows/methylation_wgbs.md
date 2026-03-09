# Workflow Pattern: Whole Genome Bisulfite Sequencing (WGBS) DNA Methylation Analysis
Category: methylation, epigenomics
Analysis type: DNA methylation, CpG methylation, bisulfite sequencing, WGBS, RRBS
Keywords: methylation, bisulfite, WGBS, RRBS, CpG, Bismark, DNA methylation, epigenomics, 5mC, differentially methylated regions, DMR

## Purpose
Quantify DNA methylation at single-CpG resolution across the genome using bisulfite
conversion followed by sequencing. Identifies differentially methylated regions (DMRs)
between conditions (e.g. cancer vs normal, treated vs control).

## Pipeline Steps (in order)
1. raw_fastqc          — FastQC on raw reads
2. trim_reads          — Trim Galore (preferred over Trimmomatic for bisulfite data):
                         adapter + quality trimming + non-CG context trimming
3. bismark_align       — Bismark: bisulfite-aware alignment to converted reference genome
4. deduplication       — bismark_deduplicate: remove PCR duplicates (critical for WGBS)
5. methylation_extract — bismark methylation extractor: extract CpG methylation calls
6. bismark2bedGraph    — Convert methylation calls to bedGraph/coverage format
7. coverage2cytosine   — Generate genome-wide CpG coverage report
8. qc_report           — Bismark2Report + MultiQC: alignment and methylation QC
9. dmr_analysis        — DSS or MethylKit (R): identify differentially methylated regions
10. dmr_annotation     — annotatr (R): annotate DMRs to CpG islands, shores, promoters

## Bisulfite-Specific Considerations
- Bismark requires a bisulfite-converted genome index (run bismark_genome_preparation once)
- Non-CG methylation: common in stem cells; use --CX flag to extract CHH/CHG contexts
- Conversion efficiency: should be >99%; check lambda spike-in or non-CpG methylation rate
- Coverage: recommend >10x per CpG for reliable methylation calls; filter low-coverage CpGs
- RRBS vs WGBS: RRBS uses restriction enzyme enrichment; trim first 2bp of R2 reads

## Key Decision Points
- Bismark vs BSMAP vs BWA-meth: Bismark is the gold standard; most widely used
- DSS vs MethylKit: DSS is more statistically rigorous; MethylKit is easier for beginners
- CpG-only vs CpHpH: for most mammalian studies, filter to CpG context only

## File Format Flow
.fastq.gz → [Trim Galore] → .fastq.gz → [Bismark] → .bam → [deduplicate] → .bam → [methylation extractor] → CpG_OT/OB.txt → [bismark2bedGraph] → .bedGraph → [DSS] → DMRs.csv

## Required Config Parameters
- samples: list of sample names
- reference: path to reference genome (Bismark will create bisulfite index)
- genome_dir: directory for Bismark genome index
- min_coverage: minimum CpG coverage to include in analysis (recommend 10)
- conditions: dict mapping sample to condition for DMR analysis

## Resource Requirements
- Bismark alignment: 32GB RAM, 8 CPUs, ~4h per sample (WGBS is slow — uses 3 CPU threads internally)
- Methylation extractor: 8GB RAM, 8 CPUs, ~1h
- DSS: 16GB RAM, 4 CPUs