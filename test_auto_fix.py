import subprocess, os
def fix_missing(pipeline_dir):
    for i in range(10):
        r = subprocess.run(["snakemake", "--dry-run", "--cores", "1"], cwd=pipeline_dir, capture_output=True, text=True)
        if r.returncode == 0: print("Passed!"); return True
        lines = (r.stdout + r.stderr).splitlines()
        missing = []
        in_affected = False
        for line in lines:
            if "affected files:" in line: in_affected=True
            elif in_affected:
                if line.startswith("        ") and line.strip(): missing.append(line.strip())
                elif line.strip(): in_affected=False
        if not missing:
            print("Failed but no affected files found:", lines[-1])
            return False
        for f in missing:
            p = os.path.join(pipeline_dir, f)
            os.makedirs(os.path.dirname(p), exist_ok=True)
            if "." in os.path.basename(f) or "json" in f or "tsv" in f or "bed" in f:
                with open(p, "w") as x: pass
            else:
                os.makedirs(p, exist_ok=True)
        print(f"Created {len(missing)} missing files/dirs. Retrying...")

