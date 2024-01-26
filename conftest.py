import pytest

from ProteinCartography import file_utils


@pytest.fixture(scope="session", autouse=True)
def repo_dirpath():
    return file_utils.find_repo_dirpath()


@pytest.fixture
def integration_test_artifacts_dirpath(repo_dirpath):
    return repo_dirpath / "ProteinCartography" / "tests" / "integration-test-artifacts"


def pytest_addoption(parser):
    """
    Add custom CLI options for pytest
    """
    parser.addoption(
        "--no-mocks",
        action="store_true",
        default=False,
        help="Run tests without mocks",
    )
