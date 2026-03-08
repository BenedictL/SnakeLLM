import json
import sys
from pathlib import Path

DEFAULT_SLURM_CONFIG = """
cluster: "sbatch --cpus-per-task={threads} --mem={resources.mem_mb} --time={resources.time}"
jobs: 50
use-conda: true
use-singularity: true
printshellcmds: true
latency-wait: 60
rerun-incomplete: true
default-resources:
  - mem_mb=4000
  - time="02:00:00"
"""

def main():
    profile_dir = Path("profiles/slurm")
    profile_dir.mkdir(parents=True, exist_ok=True)

    config_path = profile_dir / "config.yaml"
    config_path.write_text(DEFAULT_SLURM_CONFIG.strip() + "\n")

    print(f"✅ wrote {config_path}")

if __name__ == "__main__":
    main()