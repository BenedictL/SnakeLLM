"""
data/workflows/workflow_patterns.py
=====================================
Generates the Plan RAG knowledge base documents.
Run once: python -m data.workflows.workflow_patterns
Writes .md files to data/workflows/ for Plan RAG indexing.
"""

from pathlib import Path

WORKFLOWS = {

"rna_seq_de": """
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
""",


"atac_seq": """
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
""",


"wgs_variant_calling": """
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
""",


"scrna_seq": """
# Workflow Pattern: Single-Cell RNA-seq Analysis
Category: scrna-seq
Analysis type: single cell, cell clustering, cell type annotation
Keywords: scRNA-seq, single cell, 10x Genomics, Cell Ranger, Seurat, Scanpy, clustering, UMAP, cell types

## Purpose
Profile gene expression at single-cell resolution. Identify distinct cell populations,
marker genes, and cell-type-specific expression patterns.

## Pipeline Steps (in order)
1. cellranger_count    — Cell Ranger: alignment + UMI counting (10x Genomics data)
   OR star_solo        — STARsolo: open-source alternative to Cell Ranger
2. quality_control     — Seurat or Scanpy: filter low-quality cells (nCount, nFeature, % mito)
3. normalization       — Normalize + log-transform count matrix
4. feature_selection   — Select highly variable genes (HVGs)
5. dimensionality_red  — PCA → UMAP/tSNE for visualization
6. clustering          — Leiden or Louvain graph clustering
7. marker_genes        — FindMarkers / rank_genes_groups: identify cluster markers
8. cell_annotation     — Annotate clusters using known markers or automated tools (SingleR, CellTypist)
9. diff_expression     — Pseudobulk DE between conditions (DESeq2 on aggregated counts)

## Key Decision Points
- Cell Ranger vs STARsolo: Cell Ranger is the 10x gold standard; STARsolo is free and nearly equivalent
- Seurat (R) vs Scanpy (Python): Seurat more established; Scanpy integrates better with Python ML ecosystem
- Clustering resolution: test multiple resolutions (0.2 to 1.5); validate with known markers
""",

}


def write_workflow_docs():
    out_dir = Path("data/workflows")
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, content in WORKFLOWS.items():
        path = out_dir / f"{name}.md"
        path.write_text(content.strip(), encoding="utf-8")
        print(f"Written: {path}")


if __name__ == "__main__":
    write_workflow_docs()
    print("Done. Workflow pattern docs ready for Plan RAG indexing.")