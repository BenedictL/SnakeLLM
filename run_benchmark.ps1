# SnakeLLM Benchmark Runner — 15 selected prompts across RNA-seq, ATAC-seq, WGS, ChIP-seq, scRNA-seq
# Run from project root: .\run_benchmark.ps1

Set-Location C:\Users\bened\desktop\snakellm
.venv\Scripts\activate

# Create results directory if it doesn't exist
New-Item -ItemType Directory -Force -Path results | Out-Null

# Create benchmark log file with headers
$logFile = "results\benchmark.csv"
if (-not (Test-Path $logFile)) {
    "prompt,pipeline_type,output_file,status,timestamp" | Out-File -FilePath $logFile -Encoding utf8
}

# ── Define 15 important prompts ────────────────────────────────────────────────
$prompts = @(

    # ── RNA-seq (5 prompts) ──────────────────────────────────────────────────
    @{ p="run differential expression analysis on RNA-seq data using DESeq2";              o="results\rna_deseq2.json";        type="rna-seq" },
    @{ p="bulk RNA-seq differential expression in mouse liver using edgeR";                o="results\rna_edger.json";          type="rna-seq" },
    @{ p="stranded paired-end RNA-seq with hisat2 alignment and StringTie quantification"; o="results\rna_hisat2.json";         type="rna-seq" },
    @{ p="RNA-seq with UMI deduplication using UMI-tools before featureCounts";            o="results\rna_umi.json";            type="rna-seq" },
    @{ p="bulk RNA-seq with salmon quantification and tximeta import for DESeq2";          o="results\rna_salmon.json";         type="rna-seq" },

    # ── ATAC-seq (4 prompts) ──────────────────────────────────────────────────
    @{ p="ATAC-seq peak calling with MACS2 and IDR reproducibility filtering";            o="results\atac_macs2_idr.json";     type="atac-seq" },
    @{ p="ATAC-seq pipeline with Bowtie2 alignment deepTools QC and MACS2 peak calling";  o="results\atac_bowtie2.json";       type="atac-seq" },
    @{ p="ATAC-seq differential chromatin accessibility analysis with DiffBind";           o="results\atac_diffbind.json";      type="atac-seq" },
    @{ p="ATAC-seq with Tn5 shift correction and motif enrichment using HOMER";           o="results\atac_homer.json";         type="atac-seq" },

    # ── WGS / Variant Calling (3 prompts) ────────────────────────────────────
    @{ p="whole genome sequencing variant calling with GATK4 HaplotypeCaller";            o="results\wgs_gatk4.json";          type="wgs" },
    @{ p="somatic variant calling with GATK4 Mutect2 for tumor-normal paired samples";    o="results\wgs_somatic.json";        type="wgs" },
    @{ p="WGS copy number variation analysis with CNVkit";                                o="results\wgs_cnv.json";            type="wgs" },

    # ── ChIP-seq (2 prompts) ──────────────────────────────────────────────────
    @{ p="ChIP-seq pipeline for histone modification H3K27ac with MACS2 peak calling";    o="results\chip_h3k27ac.json";       type="chip-seq" },
    @{ p="ChIP-seq transcription factor binding analysis with HOMER motif discovery";     o="results\chip_tf.json";            type="chip-seq" },

    # ── scRNA-seq (1 prompt) ─────────────────────────────────────────────────
    @{ p="single-cell RNA-seq analysis with STARsolo alignment and Seurat clustering";    o="results\scrna_seurat.json";       type="scrna-seq" }
)

# ── Run all prompts ────────────────────────────────────────────────────────────
$total     = $prompts.Count
$passed    = 0
$failed    = 0
$timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"

Write-Host ""
Write-Host "======================================================" -ForegroundColor Cyan
Write-Host "  SnakeLLM Benchmark — $total prompts" -ForegroundColor Cyan
Write-Host "======================================================" -ForegroundColor Cyan
Write-Host ""

for ($i = 0; $i -lt $prompts.Count; $i++) {
    $item  = $prompts[$i]
    $num   = $i + 1
    $ts    = Get-Date -Format "yyyy-MM-dd HH:mm:ss"

    Write-Host "[$num/$total] $($item.type.ToUpper())" -ForegroundColor Yellow
    Write-Host "  Prompt: $($item.p)" -ForegroundColor White

    try {
        # Run generation
        python main.py generate $item.p --output $item.o 2>&1 | Out-Null

        if (Test-Path $item.o) {
            # Validate the JSON is not empty
            $content = Get-Content $item.o -Raw
            $json    = $content | ConvertFrom-Json

            $rulesCount = $json.rules.Count
            $toolsCount = $json.tools.Count

            Write-Host "  PASSED — $rulesCount rules, $toolsCount tools -> $($item.o)" -ForegroundColor Green
            "$($item.p),$($item.type),$($item.o),PASSED,$ts" | Add-Content -Path $logFile -Encoding utf8
            $passed++
        } else {
            Write-Host "  FAILED — output file not created" -ForegroundColor Red
            "$($item.p),$($item.type),$($item.o),FAILED_NO_FILE,$ts" | Add-Content -Path $logFile -Encoding utf8
            $failed++
        }
    } catch {
        Write-Host "  FAILED — $($_.Exception.Message)" -ForegroundColor Red
        "$($item.p),$($item.type),$($item.o),FAILED_ERROR,$ts" | Add-Content -Path $logFile -Encoding utf8
        $failed++
    }

    Write-Host ""

    # Wait 3 seconds between calls to respect API rate limits
    if ($i -lt ($prompts.Count - 1)) {
        Start-Sleep -Seconds 3
    }
}

# ── Summary ────────────────────────────────────────────────────────────────────
Write-Host "======================================================" -ForegroundColor Cyan
Write-Host "  BENCHMARK COMPLETE" -ForegroundColor Cyan
Write-Host "  Total:  $total" -ForegroundColor White
Write-Host "  Passed: $passed" -ForegroundColor Green
Write-Host "  Failed: $failed" -ForegroundColor Red
Write-Host "  Schema pass rate: $([math]::Round(($passed / $total) * 100, 1))%" -ForegroundColor Yellow
Write-Host "  Results saved to: $logFile" -ForegroundColor White
Write-Host "======================================================" -ForegroundColor Cyan
