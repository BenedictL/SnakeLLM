# SnakeLLM Benchmarking Framework

## Structure

```
benchmark/
  prompts.py    — 18 prompts (6 pipeline types × 3 difficulty levels)
  metrics.py    — metric evaluation logic
  runner.py     — orchestrates all 360 runs, caching, CSV output
  analyse.py    — summarise results into publication tables
results/
  benchmark.csv          — one row per run
  runs/
    <run_id>.json        — cached LLM response + spec per run
```

## 360 Run Matrix

| Axis          | Values                                                  | Count |
|---------------|---------------------------------------------------------|-------|
| Models        | claude-sonnet-4-6, gpt-4o, gemini-1.5-pro, qwen2.5-coder-7b | 4 |
| RAG conditions| no_rag, plan_rag_only, execute_rag_only, full_dual_rag  | 4     |
| Prompts       | 6 pipeline types × 3 difficulty levels                  | 18    |
| **Total**     |                                                         | **288** |

> Note: 4 × 4 × 18 = 288 runs. To reach 360, add a 5th model or extra prompts.

## Metrics

| Metric             | Type    | Description                                       |
|--------------------|---------|---------------------------------------------------|
| `valid_spec`       | 0/1     | Pydantic validation passes                        |
| `dry_run_pass`     | 0/1     | `snakemake --dry-run` exits 0                     |
| `rule_count_score` | 0.0–1.0 | actual / expected rule count (capped at 1.0)      |
| `tool_name_score`  | 0.0–1.0 | fraction of expected tools present (fuzzy match)  |
| `container_score`  | 0.0–1.0 | fraction of expected container substrings present |
| `dag_score`        | 0.0–1.0 | fraction of expected DAG edges present            |
| `composite_score`  | 0.0–1.0 | unweighted mean of all applicable metrics         |

## Usage

```bash
# Full benchmark (all 288/360 runs)
python -m benchmark.runner

# Single model only
python -m benchmark.runner --model claude-sonnet-4-6

# Simple prompts only (fast validation)
python -m benchmark.runner --difficulty simple

# Dry-run: print what would run, no API calls
python -m benchmark.runner --dry-run

# Re-evaluate metrics from cache (no API calls)
python -m benchmark.runner --eval-only

# Analyse results
python -m benchmark.analyse
python -m benchmark.analyse --by rag_condition
python -m benchmark.analyse --by difficulty
```

## Resuming Interrupted Runs

The runner is fully resumable. Each run's LLM response is cached to
`results/runs/<run_id>.json` immediately after the API call. If the run
is interrupted, re-running skips any run_id already present in `benchmark.csv`.

## Adding the 5th Model (fine-tuned Qwen)

Once your fine-tuned Qwen2.5-Coder-7B is ready, add it to `MODELS` in
`runner.py` and ensure `llm/providers.py` handles `"qwen2.5-coder-7b"`.
All 18 prompts × 4 RAG conditions will run automatically.
