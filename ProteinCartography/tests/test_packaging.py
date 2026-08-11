"""
Tests that the package can be installed and that its modules can be imported
outside of a git working tree.

These are unit tests and do not run the pipeline, so unlike the other tests in this
directory they do not need network access, conda environments, or mocked API responses.
"""

import ast
import os
import pathlib
import shutil
import subprocess
import sys

PACKAGE_DIRPATH = pathlib.Path(__file__).parent.parent
REPO_DIRPATH = PACKAGE_DIRPATH.parent
SETUP_FILEPATH = REPO_DIRPATH / "setup.py"


def get_setup_scripts() -> list:
    """
    Parse setup.py and return the list of filepaths in the `scripts` argument of `setup()`.

    Returns:
        The script filepaths, relative to the root of the repo.
    """
    tree = ast.parse(SETUP_FILEPATH.read_text())
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "id", None) == "setup"):
            continue
        for keyword in node.keywords:
            if keyword.arg == "scripts":
                return [element.value for element in keyword.value.elts]
    raise AssertionError("Could not find a `scripts` argument in setup.py.")


def get_module_scope_imports(filepath: pathlib.Path) -> list:
    """
    Return the names of the top-level packages imported at module scope in a module.

    Args:
        filepath: the module to parse.

    Returns:
        The top-level package names, for example "os" for `import os.path`.
    """
    module_names = []
    for node in ast.parse(filepath.read_text()).body:
        if isinstance(node, ast.ImportFrom):
            module_names.append(node.module)
        elif isinstance(node, ast.Import):
            module_names.extend(alias.name for alias in node.names)
    return [(name or "").split(".")[0] for name in module_names]


def test_all_scripts_listed_in_setup_py_exist():
    """
    Tests that every script in `scripts` exists, because setuptools raises
    a `FileNotFoundError` in its `build_scripts` step for scripts that do not,
    which makes the package impossible to install.
    """
    missing_filepaths = [
        filepath for filepath in get_setup_scripts() if not (REPO_DIRPATH / filepath).exists()
    ]
    assert not missing_filepaths, (
        f"setup.py lists scripts that do not exist: {missing_filepaths}. "
        "This makes `pip install .` fail during the `build_scripts` step."
    )


def test_no_module_imports_the_tests_package_at_module_scope():
    """
    Tests that no module in the package imports the `tests` package at module scope.

    The `tests` package is not included in the `packages` argument of `setup()`, and importing
    `tests.mocks` calls `find_repo_dirpath()` at import time, so importing it at module scope
    makes the importing module unusable outside of a git working tree.
    """
    offending_filenames = [
        filepath.name
        for filepath in sorted(PACKAGE_DIRPATH.glob("*.py"))
        if "tests" in get_module_scope_imports(filepath)
    ]
    assert not offending_filenames, (
        f"these modules import the `tests` package at module scope: {offending_filenames}. "
        "It should be imported inside the branch that uses it."
    )


def test_api_utils_is_importable_outside_a_git_repo(tmp_path):
    """
    Tests that `api_utils` can be imported from a copy of the package that is not
    inside a git working tree, which is the case for an installed copy of the package
    or an extracted release archive.
    """
    package_copy_dirpath = tmp_path / PACKAGE_DIRPATH.name
    shutil.copytree(
        PACKAGE_DIRPATH,
        package_copy_dirpath,
        ignore=shutil.ignore_patterns("integration-test-artifacts", "__pycache__"),
    )

    # `find_repo_dirpath` walks upwards, so the copy is only a valid test case if there is
    # no git working tree above it either.
    assert not any(
        (dirpath / ".git").exists()
        for dirpath in [package_copy_dirpath, *package_copy_dirpath.parents]
    )

    # `bioservices` is imported by `api_utils` but is not in the environment the tests run in,
    # so it is stubbed out; this test is about which modules `api_utils` imports, not about
    # what they contain.
    (package_copy_dirpath / "bioservices.py").write_text("class UniProt:\n    pass\n")

    # The mock-related env variables are unset because the mocks intentionally require
    # the repo, so they are not expected to work outside of a git working tree.
    env = {
        key: value for key, value in os.environ.items() if not key.startswith("PROTEINCARTOGRAPHY_")
    }
    process = subprocess.run(
        [sys.executable, "-c", "import api_utils"],
        cwd=package_copy_dirpath,
        env=env,
        capture_output=True,
        text=True,
    )
    assert process.returncode == 0, process.stderr
