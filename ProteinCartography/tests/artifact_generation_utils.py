import re

from file_utils import find_repo_dirpath
from requests.adapters import HTTPAdapter


class HTTPAdapterWithLogging(HTTPAdapter):
    """
    This class adds (crude and incomplete) logging to the `send` method

    Its intended use is to log the responses to all of the web API calls made by the pipeline
    in order to create realistic mocked API responses for use during testing.

    It is *not* intended for use in production.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._log_dirpath = find_repo_dirpath() / "logs"
        self._log_dirpath.mkdir(exist_ok=True)

    def send(self, request, **kwargs):
        response = super().send(request, **kwargs)
        self._log_response(request, response, stream=kwargs.get("stream"))
        return response

    def _log_response(self, request, response, stream=False):
        """
        log the response to a file in the self._log_dirpath directory
        """
        sanitized_url = re.sub(r"^https?://", "", request.url)
        sanitized_url = re.sub(r'[<>:"/\\|?&*]', "_", sanitized_url)
        output_filepath = self._log_dirpath / f"{request.method}__{sanitized_url[:100]}"

        if stream:
            # Consume the streamed content and write it to a file.
            content = b"".join(response.iter_content(chunk_size=128))
            with open(output_filepath, "wb") as file:
                file.write(content)

            # This allows the response content to be streamed again.
            response._content = content

        else:
            content_type = response.headers.get("Content-Type", "").split(";")[0]
            if content_type in ["application/json", "text/html", "text/plain"]:
                with open(output_filepath, "w", encoding="utf-8") as file:
                    file.write(response.text)
            else:
                with open(output_filepath, "wb") as file:
                    file.write(response.content)

        print(f"Logged response to {output_filepath}")
