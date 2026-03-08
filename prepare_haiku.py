import io

with io.open('run_benchmark.ps1', 'r', encoding='utf-8') as f:
    text = f.read()

text = text.replace('results\\benchmark.csv', 'results\\benchmark_haiku.csv')
text = text.replace('results\\', 'results\\haiku_')
# The directory is 'results' so replacing results\\ might break the directory creation
# Let's be more specific
text = text.replace('o="results\\', 'o="results\\haiku_')
text = text.replace('python main.py generate $item.p --output $item.o', 'python main.py generate $item.p --model claude-haiku-4-5 --output $item.o')
text = text.replace('SnakeLLM Benchmark -', 'SnakeLLM Benchmark (Haiku) -')

with io.open('run_benchmark_haiku.ps1', 'w', encoding='utf-8') as f:
    f.write(text)
