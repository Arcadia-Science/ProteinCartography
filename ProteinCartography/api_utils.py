#!/usr/bin/env python
import os

from bioservices import UniProt
from requests import Session
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from tests import artifact_generation_utils, mocks

__all__ = ["session_with_retry", "DefaultExpBackoffRetry", "UniProtWithExpBackoff"]

USER_AGENT_HEADER = {"User-Agent": "ProteinCartography/0.4 (Arcadia Science) python-requests/2.0.1"}

# If necessary, mock the web API responses returned by the `request` method of `requests.Session`
# Note: the env variable below is set by the `set_env_variables` pytest fixture during test setup;
# it should be used only during testing and *not* in production.
if os.environ.get("PROTEINCARTOGRAPHY_SHOULD_USE_MOCKS") == "true":
    mocks.mock_requests_session_request()

# If necessary, use the custom (and crude) `HTTPAdapterWithLogging` class in place of `HTTPAdapter`
# to log the API requests made by the pipeline.
# Note: `HTTPAdapterWithLogging` is intended for use only during development
# to help create the mocked API responses used in the pytest tests.
# The env var below should only ever be set manually by the developer in a dev environment;
# it should *not* be set in production.
if os.environ.get("PROTEINCARTOGRAPHY_SHOULD_LOG_API_REQUESTS") == "true":
    HTTPAdapter = artifact_generation_utils.HTTPAdapterWithLogging  # noqa F811


def session_with_retry():
    """
    This will return a requests Session with `DefaultExpBackoffRetry` set for the retry strategy,
    giving back a session with exponential backoff.
    """
    session = Session()
    retry = DefaultExpBackoffRetry()
    session.headers.update(USER_AGENT_HEADER)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


class DefaultExpBackoffRetry(Retry):
    """
    This class extends urllib3's `Retry` class. It sets defaults that seem to work well for the API
    calls used by the pipeline. This will give requests with the following waits: [4, 8, 16, 32, 64]
    seconds between attempts of a failing request.

    Note that if you provide an argument by position and there is a kwds value,
    either by user of default, an error will be raised by Python.
    """

    def __init__(self, *args, **kwds):
        defaults = {
            "total": 5,
            "read": 5,
            "connect": 5,
            "backoff_factor": 2,
            "status_forcelist": (500, 502, 503, 504),
            "allowed_methods": frozenset({"GET", "PUT", "POST"}),
        }

        # Update kwds to default values if not already provided.
        for key, value in defaults.items():
            kwds.setdefault(key, value)

        super().__init__(*args, **kwds)


class UniProtWithExpBackoff(UniProt):
    """
    Specialize the UniProt class to set the MAX_RETRIES to be a DefaultExpBackoffRetry
    object.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.services.settings.MAX_RETRIES = DefaultExpBackoffRetry()
