.PHONY: lint
lint:
	ruff --exit-zero check .
	snakefmt --check .

.PHONY: format
format:
	ruff --fix .
	ruff format .
	snakefmt .

.PHONY: pre-commit
pre-commit:
	pre-commit run --all-files

.PHONY: test
test:
	pytest -v

.PHONY: run-demo-workflow
run-demo-workflow:
	snakemake --snakefile Snakefile --configfile demo/config_actin.yml --use-conda $(ARGS)
