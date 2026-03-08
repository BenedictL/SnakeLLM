# enrich_benchmark.py — run from project root
import json, csv, os

input_file  = "benchmark.csv"
output_file = "benchmark_enriched.csv"

rows = []
with open(input_file, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Extract output file path from the row
        # Your CSV has it in column index 2 (pipeline_type column has the path)
        json_path = list(row.values())[2].replace("\\", "/")

        if json_path and os.path.exists(json_path):
            with open(json_path) as jf:
                spec = json.load(jf)

            rules = len(spec.get("rules", []))
            tools = len(spec.get("tools", []))

            # Check featureCounts container
            fc_ok = "N/A"
            for tool in spec.get("tools", []):
                if tool["name"].lower() in ["featurecounts", "featurecounts"]:
                    uri = tool.get("container", {}).get("full_uri", "")
                    fc_ok = "OK" if "subread" in uri else f"WRONG: {uri}"

            # Check DESeq2 aggregate rule
            deseq2_ok = "N/A"
            for rule in spec.get("rules", []):
                if "deseq2" in rule["name"].lower() or "differential" in rule["name"].lower():
                    inputs = rule.get("input", [])
                    has_wildcard = any("{sample}" in i for i in inputs)
                    deseq2_ok = "WRONG_WILDCARD" if has_wildcard else "OK"

            row["rules"]          = rules
            row["tools"]          = tools
            row["featurecounts_ok"] = fc_ok
            row["deseq2_ok"]      = deseq2_ok
        else:
            row["rules"]          = "file_not_found"
            row["tools"]          = "file_not_found"
            row["featurecounts_ok"] = "file_not_found"
            row["deseq2_ok"]      = "file_not_found"

        rows.append(row)

# Write enriched CSV
fieldnames = ["prompt","pipeline_type","model","schema_pass","rules","tools","featurecounts_ok","deseq2_ok","notes"]
with open(output_file, "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)

print(f"Saved to {output_file}")