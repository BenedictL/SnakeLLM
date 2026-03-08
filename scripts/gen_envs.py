import json
import re
import sys
from pathlib import Path

CHANNELS = ["conda-forge", "bioconda", "defaults"]

def slug(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^a-z0-9._-]+", "-", s)
    s = re.sub(r"-+", "-", s).strip("-")
    return s or "env"

def env_name_from_spec(spec_path: Path) -> str:
    return slug(spec_path.stem)

def tool_to_dep(tool: dict) -> str | None:
    name = tool.get("name")
    ver = tool.get("version")
    if not name:
        return None
    name = name.strip()
    if ver and str(ver).strip() and str(ver).strip().lower() not in {"latest", "*"}:
        return f"{name}={ver}"
    return name

def main():
    if len(sys.argv) < 2:
        print("Usage: python gen_envs.py specs/*.json")
        sys.exit(1)

    out_dir = Path("envs")
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in sys.argv[1:]:
        spec_path = Path(f)
        data = json.loads(spec_path.read_text())

        tools = data.get("tools", [])
        deps = []
        for t in tools:
            dep = tool_to_dep(t)
            if dep:
                deps.append(dep)

        # de-dup while preserving order
        seen = set()
        deps = [d for d in deps if not (d in seen or seen.add(d))]

        env_name = env_name_from_spec(spec_path)
        env_file = out_dir / f"{env_name}.yaml"

        yaml = []
        yaml.append(f"name: {env_name}")
        yaml.append("channels:")
        for c in CHANNELS:
            yaml.append(f"  - {c}")
        yaml.append("dependencies:")
        yaml.append("  - python=3.11")
        for d in deps:
            yaml.append(f"  - {d}")
        yaml.append("")

        env_file.write_text("\n".join(yaml))
        print(f"✅ wrote {env_file} ({len(deps)} tool deps)")

if __name__ == "__main__":
    main()