.PHONY: lint
lint:
	ruff check --exit-zero .
	snakefmt --check .

.PHONY: format
format:
	ruff format .
	ruff check --fix .
	snakefmt .

.PHONY: pre-commit
pre-commit:
	pre-commit run --all-files

.PHONY: test
test:
	pytest -vv -s .

.PHONY: test-without-mocks
test-without-mocks:
	pytest -vv -s --no-mocks

.PHONY: run-demo-workflow
run-demo-workflow:
	snakemake \
		--snakefile Snakefile \
		--configfile demo/search-mode/config_actin.yml \
		--use-conda  \
		--cores all \
		$(ARGS)

.PHONY: generate-search-mode-rulegraph
generate-search-mode-rulegraph:
	snakemake \
		--configfile demo/search-mode/config_actin.yml \
		--rulegraph \
		| dot -Tpng \
		> rulegraph-search-mode.png

.PHONY: generate-cluster-mode-rulegraph
generate-cluster-mode-rulegraph:
	snakemake \
		--configfile demo/cluster-mode/config.yml \
		--rulegraph \
		| dot -Tpng \
		> rulegraph-cluster-mode.png
