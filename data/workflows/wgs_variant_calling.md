# Workflow Pattern: Whole Genome Sequencing Variant Calling (GATK Best Practices)
Category: wgs, variant-calling
Analysis type: SNP calling, indel calling, germline variants, somatic variants
Keywords: WGS, variant calling, GATK, SNP, indel, germline, somatic, BWA, HaplotypeCaller, VCF, genotyping

## Purpose
Identify genetic variants (SNPs, indels, structural variants) from whole genome or whole exome sequencing data.
Follows GATK4 Best Practices pipeline for germline short variant discovery.

## Pipeline Steps (in order)
1. raw_fastqc            — FastQC on raw reads
2. trim_reads            — Trimmomatic: adapter and quality trimming
3. align_bwa             — BWA-MEM2: align reads to reference genome (faster than BWA-MEM)
4. sort_bam              — samtools sort: coordinate sort aligned BAM
5. mark_duplicates       — Picard MarkDuplicates: flag PCR duplicates (do NOT remove for GATK)
6. base_recalibration    — GATK BaseRecalibrator + ApplyBQSR: correct systematic sequencing errors
7. haplotype_caller      — GATK HaplotypeCaller: call variants per sample in GVCF mode (-ERC GVCF)
8. combine_gvcfs         — GATK CombineGVCFs or GenomicsDBImport: merge per-sample GVCFs
9. genotype_gvcfs        — GATK GenotypeGVCFs: joint genotyping across all samples
10. variant_filtration   — GATK VQSR (>30 samples) or hard filtering (<30 samples)
11. annotate_variants    — GATK Funcotator or VEP: functional annotation of variants
12. variant_qc           — bcftools stats + MultiQC: summary statistics on VCF

## Key Decision Points
- BWA vs BWA-MEM2: BWA-MEM2 is 2-3x faster and drop-in compatible; always prefer it
- VQSR vs hard filtering: VQSR requires >30 WGS samples or >10 WES samples; use hard filters otherwise
- WGS vs WES: add interval_list parameter for WES to restrict to capture regions
- Germline vs somatic: this pattern is germline; for somatic (tumor/normal) use GATK Mutect2 instead

## GATK-Specific Requirements
- Reference genome must be indexed: samtools faidx + picard CreateSequenceDictionary
- Known sites VCFs required for BQSR: dbSNP, 1000G indels, Mills indels
- For joint calling: ALL samples must be called with HaplotypeCaller -ERC GVCF first

## File Format Flow
.fastq.gz → [BWA-MEM2] → .bam → [sort+markdup+BQSR] → recal.bam → [HaplotypeCaller] → .g.vcf.gz → [joint genotyping] → .vcf.gz → [filter+annotate] → annotated.vcf.gz

## Required Config Parameters
- samples: list of sample names
- reference: path to reference genome (must be indexed)
- known_sites: list of paths to known variant VCFs (dbSNP, 1000G)
- intervals: (optional) BED file for WES target capture regions
- scatter_count: number of scatter intervals for parallelization (recommend 24 for WGS)

## Resource Requirements (intensive pipeline)
- BWA-MEM2: 32GB RAM, 16 CPUs, ~3h per 30x WGS sample
- HaplotypeCaller: 16GB RAM, 4 CPUs, ~8h per 30x WGS sample (use scatter-gather)
- GenotypeGVCFs: 32GB RAM, 8 CPUs