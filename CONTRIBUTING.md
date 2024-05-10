# Contributing to ProteinCartography
Thanks for your interest in contributing to ProteinCartography!
Please read this document in its entirety before contributing to help ensure that your contribution meets our standards and is readily accepted.

## Getting Started
All the packages needed to develop for ProteinCartography are found in the `envs/cartography_dev.yml` conda environment.
You can install this environment as follows:

1. Make sure `miniconda` is installed. Even if you’re using an Apple Silicon (M1, M2, etc. macOS) laptop, you will need to install the macOS Intel x86-64 version of `miniconda` [here](https://docs.conda.io/projects/miniconda/en/latest/).

2. Create a conda environment from the `cartography_dev.yml` file in the `envs/` directory.
```sh
conda env create -n cartography_dev --file envs/cartography_dev.yml
```

3. Activate the environment.
```sh
conda activate cartography_dev
```

## How to contribute
### Bug reports and feature requests
We track all bugs, new feature requests, enhancements, etc. using [GitHub Issues](https://github.com/Arcadia-Science/ProteinCartography/issues). Please check to make sure that your issue has not already been reported. If it has, please add a comment to the existing issue instead of creating a new one.

### Making a contribution
The steps below apply to both external and internal contributors and also apply to working both with this repo itself and with your own fork. However, if you are an external contributor, please fork this repository and make your changes in a new branch on your fork. Please see the GitHub docs on [working with forks](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) for more details about how to do this.

1. Whenever you start work, you should make sure to `git pull` on the main branch to get the latest version of the `main` branch.

1. When you’re planning on working on an issue, you should claim it by assigning it to yourself. You should also add a comment briefly explaining how you plan to address the issue. If you are an external contributor, please wait for a maintainer to sign off on your plan before you start working; this will make it easier for us to accept your PR later on.

1. After claiming an issue, create a new branch from the `main` branch. __Make sure your branch begins with your initials, followed by a forward slash__. This is very important to keep everyone's branches well-organized. Use short descriptive names for the branch name itself. For example, if your initials are `abc` and you are adding a new feature to evaluate clustering, you might name your branch `abc/add-cluster-evaluation`. If you are working on an issue to fix a bug, you might name your branch `abc/fix-foldseek-format-bug`.
	Create the new branch using the following command:
	```sh
	git checkout -b <your-initials>/<branch-name>
	```

1. Once you’ve created a branch, push the branch to the GitHub repo so you can keep a remote record of your work.
	```sh
	git push -u origin <your-initials>/<branch-name>
	```

1. Once you’ve completed the feature or fixed the bug, you are ready to open a PR. Please aim to keep the PRs as small as possible to increase the speed and ease of review; a PR can be as small as changing a few characters or resolving a single bug. When you open a PR, please use a succinct and human-readable title and always add a description of the changes you made as well as a link to the issue that the PR addresses.

1. Check that your PR is passing the CI checks. These checks verify that the changes in your PR are formatted correctly (we use a tool called `ruff` for formatting; see below for details) and also that your PR passes the automated tests of the pipeline.

### Keeping your development branches up to date
Occasionally, your development branch will be behind the main branch. If this happens and you’re working on a file that was changed in the updated version of the main branch, you may need to merge the updated `main` branch into your local development branch.

1. First, update the main branch in your local repo using the following in main:
	```sh
	git checkout main
	git pull
	```

1. Next, check out your development branch:
	```sh
	git checkout <your-initials>/<branch-name>
	```

1. Now, merge the main branch into your local branch:
	```sh
	git merge origin/main
	```
	__Note that you must be on your local development branch when you call `git merge`!__

1. Once you’ve merged the main branch into your local development branch, use `git push` to push the merged changes to your branch on GitHub.

### Linting
We use `ruff` to lint and format our Python code. We also use snakemake’s code formatter, `snakefmt` , to format the snakefiles. You can and should run these tools in your local repo using the commands `make format` and `make lint`. Note that `ruff` is also available as an extension in VS Code, allowing you to configure VS Code to automatically format your code whenever you save the file you are currently editing.

### Testing
Tests are found in the `ProteinCartography/tests/` directory. We use `pytest` for testing; you can run the tests locally using the command `make test`. Currently, we only have integration-like tests that run the pipeline in both 'search' and 'cluster' modes using a test dataset and test config designed to allow the pipeline run very quickly (2-3min). The tests then check that the output files are created and have the correct shape. We plan to add unit tests in the future.

#### Running the tests without mocked API responses
When the pipeline is run in 'search' mode, it makes many calls to external APIs (including Foldseek, Blast, and Alphafold). By default, these calls are mocked during testing so that the tests do not depend on external APIs; this is important to ensure test reproducibility and also helps to make the tests run quickly. However, it is important to periodically test that the pipeline also runs correctly when real API calls are made. To do this, you can run the tests without mocks using `make test-without-mocks`.

When merging PRs on GitHub, it is likewise important to test that the pipeline runs correctly with real API calls. To do so, add the label `run-slow-tests` to your PR. This will trigger the CI actions (see below) to run again on your PR, but now without mocks. __Please add this label only when your PR is ready to merge, as it will cause the CI to run more slowly and will also result in unnecessary API calls.__

#### Updating the mocked API responses
When changes you have made involve changes to the API calls made by the pipeline, it will be necessary to update the mocked responses in order for the tests to pass. Currently, this is a manual process.
1. Enable API response logging in your local environment by setting the following environment variable:
	```sh
	export PROTEINCARTOGRAPHY_SHOULD_LOG_API_REQUESTS=true
	```

1. Run the pipeline in 'search' mode using the 'small' search-mode demo (this demo uses the same input PDB file as the tests). The API responses made by the pipeline will be logged a `logs/` directory in the root of the repo.
	```sh
	snakemake --cores all --use-conda --configfile demo/search-mode/config_actin_small.yml
	```

1. Use the logged responses to update the mocked responses constructed in the `ProteinCartography/tests/mocks.py` module. For large responses, response contents should be added to `ProteinCartography/tests/integration-test-artifacts/search-mode/actin/api_response_content/`.

1. When you're finished, don't forget to delete the `logs/` directory and unset the `PROTEINCARTOGRAPHY_SHOULD_LOG_API_REQUESTS` environment variable.

### CI
We use GitHub Actions for CI. Currently, there is one workflow for linting and one for testing. Both workflows are run automatically on GitHub when a PR is opened and also whenever new commits are pushed to an open PR. PRs cannot be merged until the CI checks pass.

The linting workflow runs `ruff --check` and `snakefmt --check` on all Python and snakefiles in the repo. This means that the workflow does not modify any files; it only checks that your code is formatted correctly. If the workflow fails for your PR, you can run `make format` locally to format your code and `make lint` to determine if there are lint errors that need to be fixed.

The testing workflow runs pytest using the same `make test` command that is used locally. If the workflow fails for your PR, it is usually best to run `make test` locally to recapitulate the failure and determine which tests are failing.

### Style guide
In addition to the formatting and lint rules imposed by `ruff` and `snakefmt`, we also have a few additional style rules that are not enforced by these tools. These rules are listed below.

- Function and variable names should be in `lower_snake_case` and should be descriptive; avoid abbreviations.
- Function arguments and return values should have [type hints](https://docs.python.org/3/library/typing.html).
- Functions should include a [Google-style docstring](https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings) explaining the arguments and what the function returns (if not `None`).
- Comments should be written in complete sentences in the present tense and should end with a period.
- Comments should be used sparingly and only when necessary to explain something that is not obvious from the code itself.
- Class names should use `CapitalizedCamelCase` with descriptive names.
- Currently, we don’t use many custom classes, but the conventions for functions apply to class methods as well.

Here is an example of a function that adheres to all of these rules:
```python
def add_integers(first_integer: int, second_integer: int) -> int:
	"""
	Add two integers together, returning the result.

	Args:
		first_integer (int): first integer to add.
		second_integer (int): second integer to add.

	Returns:
		The sum of the two integers.
	"""
	result = first_integer + second_integer
	return result
```

### Code organization
We strive to encapsulate new functionality within modular Python scripts that accept arguments from the command line using `argparse`. These scripts are then called from snakemake rules and can also be run directly from the command line by the user.
- Every script should include a `parse_args()` function and a `main()` function.
- Every script with `#!/usr/bin/env python` (so that the scripts are executable from the command line on unix systems).
- An example template for new scripts is found in [`template.py`](./ProteinCartography/template.py).

### Adding new dependencies
First, please consider carefully whether you need to add a new dependency to the project.
When changes you have made absolutely require new dependencies, please make sure that they are `conda`-installable.
Dependencies should be added to two environment files:
1. the `cartography_dev.yml` file in the `envs/` directory.
2. the appropriate snakemake rule environment file in the `envs/` directory.

In both files, please include the version of the dependency you are using (this is called "pinning" the dependency).
Include only the exact version number; do not include the package hash.
For example, if you are adding a new dependency called `new_dependency` and you are using version `1.2.3`, you would add the following line to the `cartography_dev.yml` file:
```yaml
- new_dependency=1.2.3
```


## Crediting Contributions
See how we recognize feedback and contributions to our code at Arcadia [here](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
