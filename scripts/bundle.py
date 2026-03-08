import json
import shutil
import sys
from pathlib import Path
from zipfile import ZipFile

def main():
    if len(sys.argv) < 2:
        print("Usage: python bundle.py specs/*.json")
        sys.exit(1)

    dist_dir = Path("dist")
    dist_dir.mkdir(exist_ok=True)

    for spec_file in sys.argv[1:]:
        spec_path = Path(spec_file)
        name = spec_path.stem

        zip_path = dist_dir / f"{name}.zip"

        with ZipFile(zip_path, "w") as z:
            # include spec
            z.write(spec_path, arcname=f"{name}/spec.json")

            # include env if exists
            env_file = Path("envs") / f"{name}.yaml"
            if env_file.exists():
                z.write(env_file, arcname=f"{name}/env.yaml")

            # include slurm profile
            profile_file = Path("profiles/slurm/config.yaml")
            if profile_file.exists():
                z.write(profile_file, arcname=f"{name}/profiles/slurm/config.yaml")

        print(f"✅ wrote {zip_path}")

if __name__ == "__main__":
    main()