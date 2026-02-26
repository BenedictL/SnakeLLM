"""
llm/biocontainers_indexer.py
============================
Fetches tool + container metadata from the BioContainers API and builds
a searchable index for the Execute RAG layer.

BioContainers API docs: https://api.biocontainers.pro/ga4gh/trs/v2/tools

Run this script ONCE to populate ./data/biocontainers/
Then Execute RAG reads from there — no live API calls during inference.

Usage:
    python -m llm.biocontainers_indexer --tools star trimmomatic samtools deseq2
    python -m llm.biocontainers_indexer --from-file data/tool_list.txt
"""

import argparse
import json
import time
import logging
from pathlib import Path
from typing import Optional
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)s  %(message)s")
log = logging.getLogger(__name__)

BIOCONTAINERS_API = "https://api.biocontainers.pro/ga4gh/trs/v2/tools"
OUTPUT_DIR        = Path("data/biocontainers")
RATE_LIMIT_SECS   = 0.5   # be polite to the API


# ── FETCH FROM API ────────────────────────────────────────────────────────────

def fetch_tool_metadata(tool_name: str) -> Optional[dict]:
    """
    Queries the BioContainers TRS v2 API for a single tool.
    Returns the raw API response dict, or None on failure.
    """
    url    = f"{BIOCONTAINERS_API}/{tool_name}"
    params = {"format": "json"}

    try:
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        return resp.json()
    except requests.HTTPError as e:
        if resp.status_code == 404:
            log.warning(f"Tool '{tool_name}' not found in BioContainers registry")
        else:
            log.error(f"HTTP error fetching '{tool_name}': {e}")
        return None
    except requests.RequestException as e:
        log.error(f"Network error fetching '{tool_name}': {e}")
        return None


def fetch_tool_versions(tool_name: str) -> list[dict]:
    """
    Fetches all available versions/containers for a tool.
    Returns list of version dicts sorted by recency.
    """
    url = f"{BIOCONTAINERS_API}/{tool_name}/versions"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        versions = resp.json()
        # Sort by version tag — newest first (rough lexicographic sort)
        return sorted(versions, key=lambda v: v.get("name", ""), reverse=True)
    except Exception as e:
        log.warning(f"Could not fetch versions for '{tool_name}': {e}")
        return []


def search_biocontainers(query: str, limit: int = 5) -> list[dict]:
    """
    Text search across BioContainers registry.
    Useful for finding containers when exact tool name is unknown.
    """
    params = {
        "name":   query,
        "limit":  limit,
        "offset": 0,
        "format": "json"
    }
    try:
        resp = requests.get(BIOCONTAINERS_API, params=params, timeout=10)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        log.error(f"Search failed for '{query}': {e}")
        return []


# ── PARSE INTO STRUCTURED FORMAT ──────────────────────────────────────────────

def parse_container_record(tool_name: str, raw: dict, versions: list[dict]) -> dict:
    """
    Converts raw BioContainers API response into a clean structured record
    that the Execute RAG system can embed and retrieve.
    """
    # Extract best container URI — prefer quay.io Singularity, fall back to Docker
    best_container = extract_best_container(tool_name, versions)

    # Extract description
    description = raw.get("description", "").strip()
    if not description and versions:
        description = versions[0].get("description", "")

    # Extract all available tags
    all_tags = []
    for v in versions[:10]:  # limit to 10 most recent
        for img in v.get("images", []):
            tag = img.get("image_name", "")
            if tag:
                all_tags.append({
                    "tag":      tag,
                    "type":     img.get("image_type", ""),  # Docker, Singularity
                    "registry": img.get("registry", ""),
                })

    return {
        "tool_name":       tool_name,
        "display_name":    raw.get("name", tool_name),
        "description":     description,
        "homepage":        raw.get("url", ""),
        "tool_class":      raw.get("tool_class", {}).get("name", ""),
        "best_container":  best_container,
        "all_tags":        all_tags,
        "version_count":   len(versions),
        # Metadata for Execute RAG embedding
        "embedding_text":  build_embedding_text(tool_name, description, best_container, versions),
    }


def extract_best_container(tool_name: str, versions: list[dict]) -> dict:
    """
    Priority order for container selection:
    1. quay.io Docker image (most common for BioContainers)
    2. docker.io image
    3. Any available image
    Returns empty dict if none found.
    """
    for v in versions[:5]:  # check 5 most recent versions
        images = v.get("images", [])

        # Prefer quay.io
        for img in images:
            name = img.get("image_name", "")
            if "quay.io" in name and img.get("image_type") == "Docker":
                return {
                    "full_uri": name,
                    "registry": "quay.io",
                    "tag":      name.split(":")[-1] if ":" in name else "latest",
                    "version":  v.get("name", ""),
                    "type":     "Docker",
                }

        # Fall back to docker.io
        for img in images:
            name = img.get("image_name", "")
            if img.get("image_type") == "Docker":
                return {
                    "full_uri": name,
                    "registry": name.split("/")[0] if "/" in name else "docker.io",
                    "tag":      name.split(":")[-1] if ":" in name else "latest",
                    "version":  v.get("name", ""),
                    "type":     "Docker",
                }

    return {}


def build_embedding_text(
    tool_name: str,
    description: str,
    best_container: dict,
    versions: list[dict]
) -> str:
    """
    Builds the text that gets embedded into the vector DB.
    Designed to be retrievable by semantic queries about the tool's purpose.
    """
    recent_versions = [v.get("name", "") for v in versions[:3]]
    container_uri   = best_container.get("full_uri", "not available")

    return f"""
Tool: {tool_name}
Description: {description}
Container: {container_uri}
Recent versions: {', '.join(recent_versions)}
Usage context: This tool is used in bioinformatics pipelines for {description.lower()}.
Snakemake container directive: container: "{container_uri}"
""".strip()


# ── SAVE TO DISK ──────────────────────────────────────────────────────────────

def save_tool_record(record: dict) -> Path:
    """Saves a parsed tool record as JSON to the output directory."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / f"{record['tool_name']}.json"
    with open(path, "w") as f:
        json.dump(record, f, indent=2)
    log.info(f"  Saved → {path}")
    return path


# ── MAIN INDEXING PIPELINE ────────────────────────────────────────────────────

def index_tool(tool_name: str) -> Optional[dict]:
    """Full pipeline: fetch → parse → save for one tool."""
    log.info(f"Indexing: {tool_name}")

    raw      = fetch_tool_metadata(tool_name)
    versions = fetch_tool_versions(tool_name)

    if raw is None and not versions:
        log.warning(f"  Skipping '{tool_name}' — no data found")
        return None

    raw      = raw or {}
    record   = parse_container_record(tool_name, raw, versions)
    save_tool_record(record)

    log.info(f"  ✓ {tool_name} → {record['best_container'].get('full_uri', 'no container found')}")
    return record


def index_all(tool_names: list[str]) -> dict[str, dict]:
    """Index a list of tools with rate limiting."""
    results = {}
    for i, name in enumerate(tool_names):
        record = index_tool(name.lower().strip())
        if record:
            results[name] = record
        if i < len(tool_names) - 1:
            time.sleep(RATE_LIMIT_SECS)
    return results


# ── DEFAULT TOOL LIST ─────────────────────────────────────────────────────────
# These are the tools covering RNA-seq, ATAC-seq, and WGS pipelines.

DEFAULT_TOOLS = [
    # QC
    "fastqc", "multiqc", "fastp",
    # Trimming
    "trimmomatic",
    # RNA-seq alignment
    "star", "hisat2",
    # WGS/WES alignment
    "bwa", "bwa-mem2",
    # BAM processing
    "samtools", "picard",
    # RNA-seq quantification
    "featurecounts", "htseq",
    # Variant calling
    "gatk4", "bcftools", "deepvariant",
    # ATAC-seq peak calling
    "macs2", "macs3",
    # ATAC-seq QC
    "deeptools", "ataqv",
    # Genome annotation
    "bedtools", "bedops",
    # R / Bioconductor (these have containers too)
    "bioconductor-deseq2", "bioconductor-edger",
    "bioconductor-clusterprofiler",
    "bioconductor-chipseeker",
    # Methylation
    "bismark",
    # Misc
    "bowtie2", "subread",
    
    # Newly added extended tools
    "salmon", "kallisto", "stringtie", "snpeff",
    "manta", "cnvkit", "rseqc", "qualimap",
    "starsolo", "cellranger", "diffbind"
]


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Index BioContainers tool metadata for Execute RAG")
    group  = parser.add_mutually_exclusive_group()
    group.add_argument("--tools",     nargs="+", help="Space-separated tool names to index")
    group.add_argument("--from-file", type=str,  help="Text file with one tool name per line")
    group.add_argument("--defaults",  action="store_true", help="Index the default tool list")
    args = parser.parse_args()

    if args.tools:
        tools = args.tools
    elif args.from_file:
        tools = Path(args.from_file).read_text().strip().splitlines()
    else:
        tools = DEFAULT_TOOLS

    log.info(f"Starting BioContainers indexing for {len(tools)} tools...")
    results = index_all(tools)
    log.info(f"\n✓ Done. Indexed {len(results)}/{len(tools)} tools → {OUTPUT_DIR}/")