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
Keywords: scRNA-seq, single cell, 10x Genomics, Cell Ranger, cellranger, Seurat, Scanpy, clustering, UMAP, cell types

## Purpose
Profile gene expression at single-cell resolution. Identify distinct cell populations,
marker genes, and cell-type-specific expression patterns.

## Pipeline Steps (in order)
1. cellranger_count    — Cell Ranger (cellranger): alignment + UMI counting (10x Genomics data)
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


# ── NEW: ChIP-seq ─────────────────────────────────────────────────────────────
# plan_rag.py has synonym expansion for "chip-seq" but no doc existed before.
# Without this, any ChIP-seq prompt returns 0 Plan RAG hits and the LLM
# falls back to ATAC-seq docs (wrong assay, wrong tools).

"chip_seq": """
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
""",


# ── NEW: WGBS Methylation ──────────────────────────────────────────────────────
# Methylation queries have synonym expansion in plan_rag.py but no doc existed.

"methylation_wgbs": """
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