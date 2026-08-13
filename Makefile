SNAKEFILE := workflow/Snakefile
LOCAL_PROFILE := profiles/local
SNAKEMAKE := uv run snakemake --snakefile $(SNAKEFILE)
SNAKEMAKE_ARGS ?=

.DEFAULT_GOAL := help

.PHONY: help install inputs validate-inputs dry-run reproduce lint test check

help: ## Show the available commands
	@awk 'BEGIN {FS = ":.*## "; printf "Usage: make <target>\n\n"} /^[a-zA-Z_-]+:.*## / {printf "  %-16s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install the exact locked environment
	uv sync --frozen

inputs: ## Download and verify the openly redistributable inputs
	./scripts/fetch_open_inputs.sh

validate-inputs: ## Validate the complete frozen input tree
	$(SNAKEMAKE) --profile $(LOCAL_PROFILE) results/provenance/input_manifest.json $(SNAKEMAKE_ARGS)

dry-run: ## Resolve the complete workflow without executing jobs
	$(SNAKEMAKE) --profile $(LOCAL_PROFILE) --dry-run --forceall $(SNAKEMAKE_ARGS)

reproduce: ## Rebuild all analyses, figures, tables, and audits
	$(SNAKEMAKE) --profile $(LOCAL_PROFILE) $(SNAKEMAKE_ARGS)

lint: ## Run static checks
	uv run --group dev ruff check src tests

test: ## Run the unit test suite
	uv run --group dev pytest -q

check: lint test ## Run all fast quality gates
	uv lock --check
