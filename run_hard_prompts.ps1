# run_hard_prompts.ps1
Set-Location C:\Users\bened\desktop\snakellm

$prompts = @(
    @{ p="paired-end RNA-seq with UMI deduplication using UMI-tools followed by STAR alignment featureCounts quantification and DESeq2 differential expression with batch correction using limma"; o="hard_results/hard_01_umi_batch.json" },
    @{ p="single-end and paired-end mixed ATAC-seq cohort with automatic detection of library type Bowtie2 alignment mitochondrial read filtering Picard duplicate removal MACS2 peak calling IDR reproducibility filtering and DiffBind differential accessibility"; o="hard_results/hard_02_mixed_atac.json" },
    @{ p="whole genome sequencing tumor-normal paired somatic variant calling with BWA-MEM2 alignment GATK4 Mutect2 somatic SNV calling Manta structural variant detection CNVkit copy number analysis and VEP variant annotation"; o="hard_results/hard_03_tumor_normal_wgs.json" },
    @{ p="single-cell RNA-seq with Cell Ranger alignment ambient RNA removal using SoupX doublet detection with scDblFinder Seurat clustering cell type annotation with SingleR and trajectory analysis with Monocle3"; o="hard_results/hard_04_scrna_full.json" },
    @{ p="paired-end ChIP-seq for three histone marks H3K4me3 H3K27ac and H3K27me3 with Bowtie2 alignment Picard duplicate marking MACS2 narrow and broad peak calling IDR filtering deepTools bigWig generation and ChIPseeker annotation"; o="hard_results/hard_05_multihiston_chip.json" },
    @{ p="bulk RNA-seq with both STAR featureCounts and salmon alevin quantification paths running in parallel followed by DESeq2 on featureCounts output and tximeta DESeq2 on salmon output with comparison of both results using MultiQC"; o="hard_results/hard_06_dual_quant.json" },
    @{ p="long-read Oxford Nanopore RNA-seq with minimap2 splice-aware alignment featureCounts quantification DESeq2 differential expression and comparison against short-read Illumina results using the same sample set"; o="hard_results/hard_07_nanopore_rnaseq.json" },
    @{ p="ATAC-seq and RNA-seq multi-omics integration on matched samples with separate preprocessing pipelines followed by chromatin accessibility peak annotation overlapping differentially expressed genes and transcription factor motif enrichment at open chromatin regions near DEGs"; o="hard_results/hard_08_multiomics_atac_rna.json" },
    @{ p="bisulfite sequencing whole genome methylation analysis with Bismark alignment deduplication methylation extraction DMR calling using DSS CpG island annotation and integration with RNA-seq gene expression data from matched samples"; o="hard_results/hard_09_wgbs_rnaseq.json" },
    @{ p="germline trio variant calling with BWA-MEM2 alignment GATK4 HaplotypeCaller per-sample GVCF calling joint genotyping across mother father and proband VQSR filtering DeNovoGear de novo variant detection and SnpEff functional annotation"; o="hard_results/hard_10_trio_wgs.json" }
)

$passed = 0
$failed = 0
$ts     = Get-Date -Format "yyyy-MM-dd HH:mm:ss"

New-Item -ItemType Directory -Force -Path "hard_results" | Out-Null
New-Item -ItemType Directory -Force -Path "hard_outputs" | Out-Null

"prompt,output_file,status,timestamp" | Out-File hard_results\hard_benchmark.csv -Encoding utf8

foreach ($item in $prompts) {
    $p_short = $item.p.Substring(0,60)
    Write-Host "Generating: $p_short..." -ForegroundColor Cyan
    try {
        python main.py generate $item.p --output $item.o 2>&1 | Out-Null
        if (Test-Path $item.o) {
            $json = Get-Content $item.o | ConvertFrom-Json
            $rc = $json.rules.Count
            $tc = $json.tools.Count
            Write-Host "  PASSED - $rc rules, $tc tools" -ForegroundColor Green
            $csv_line = "{0},{1},PASSED,{2}" -f $item.p, $item.o, $ts
            $csv_line | Add-Content hard_results\hard_benchmark.csv -Encoding utf8
            $passed++

            Write-Host "  Building architecture in hard_outputs..." -ForegroundColor Cyan
            python generate_snakefile.py $item.o hard_outputs
        } else {
            Write-Host "  FAILED - no output file" -ForegroundColor Red
            $csv_line = "{0},{1},FAILED,{2}" -f $item.p, $item.o, $ts
            $csv_line | Add-Content hard_results\hard_benchmark.csv -Encoding utf8
            $failed++
        }
    } catch {
        Write-Host "  ERROR - $($_.Exception.Message)" -ForegroundColor Red
        $failed++
    }
    Start-Sleep -Seconds 3
}

Write-Host ""
$total = $passed + $failed
Write-Host "Hard prompt results: $passed passed / $total total" -ForegroundColor Yellow
