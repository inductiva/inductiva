"""Methods to interact with the available simulators."""
import requests
from typing import Any
import logging


# Constants
DEFAULT_API_URL: str = "https://api-dev.inductiva.ai/simulators/available-images"
DEFAULT_TIMEOUT: int = 120

ERROR_HTTP: str = "API returned HTTP error: {}"
ERROR_CONNECTION: str = "Failed to connect to the API."
ERROR_TIMEOUT: str = "The request timed out."
ERROR_REQUEST: str = "An unexpected error occurred while making the API request."


def list_available_images() -> Any:
    """Fetch the list of available simulator images from the API."""
    try:
        response = requests.get(DEFAULT_API_URL, timeout=DEFAULT_TIMEOUT)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"{ERROR_HTTP} {http_err} - Status code: {response.status_code}")
        raise RuntimeError(ERROR_HTTP.format(http_err))
    except requests.exceptions.ConnectionError  as conn_err:
        logging.error(f"{ERROR_CONNECTION}: {conn_err}")
        raise RuntimeError(ERROR_CONNECTION)
    except requests.exceptions.Timeout as timeout_err:
        logging.error(f"{ERROR_TIMEOUT}: {timeout_err}")
        raise RuntimeError(ERROR_TIMEOUT)
    except requests.exceptions.RequestException as req_err:
        logging.error(f"{ERROR_REQUEST}: {req_err}")
        raise RuntimeError(ERROR_REQUEST)

