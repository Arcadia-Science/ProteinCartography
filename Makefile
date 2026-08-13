# These checks are the same ones that `.github/workflows/lint.yml` runs, in the same order,
# so that a passing `make lint` means a passing lint workflow.
.PHONY: lint
lint:
	ruff check .
	ruff format --check .
	snakefmt --check .

# Note that `ruff check --fix` runs before `ruff format`, because the fixer's edits are not
# themselves formatted; in the other order, a fix can leave the tree failing `make lint`.
.PHONY: format
format:
	ruff check --fix .
	ruff format .
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
