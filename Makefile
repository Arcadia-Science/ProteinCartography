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

.PHONY: test-without-mocks
test-without-mocks:
	pytest -vv -s --no-mocks

.PHONY: run-demo-workflow
run-demo-workflow:
	snakemake --snakefile Snakefile --configfile demo/config_actin.yml --use-conda $(ARGS)
