PY ?= python3
SPECS := $(wildcard specs/*.json)

.PHONY: help validate envs profile bundle all clean

help:
	@echo "Targets:"
	@echo "  make validate   # validate all specs"
	@echo "  make envs       # generate env yamls"
	@echo "  make profile    # generate slurm profile"
	@echo "  make bundle     # bundle zip artifacts"
	@echo "  make all        # validate + envs + profile + bundle"
	@echo "  make clean      # remove generated artifacts"

validate:
	$(PY) scripts/validate_spec.py $(SPECS)

envs:
	$(PY) scripts/gen_envs.py $(SPECS)

profile:
	$(PY) scripts/gen_slurm_profile.py

bundle:
	$(PY) scripts/bundle.py $(SPECS)

all: validate envs profile bundle

clean:
	rm -rf dist/*.zip