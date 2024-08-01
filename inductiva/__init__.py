"""Client for Inductiva's web API."""
import os
import sys
import logging
import contextvars
from urllib3.exceptions import MaxRetryError, NewConnectionError

from inductiva.client.apis.tags.version_api import VersionApi
from inductiva.client.configuration import Configuration
from inductiva.client.exceptions import ApiException
from inductiva._cli.cmd_user.info import get_info
from inductiva.client.api_client import ApiClient

from . import simulators
from . import resources
from . import projects
from . import commands
from . import storage
from . import utils
from . import tasks
from . import users
from . import logs
from . import api
from .templating import TemplateManager

logs.setup(getattr(logging, os.environ.get("INDUCTIVA_LOG_LEVEL", "INFO")))

api_url = os.environ.get("INDUCTIVA_API_URL", "https://api.inductiva.ai")
_output_dir = contextvars.ContextVar("INDUCTIVA_OUTPUT_DIR",
                                     default=os.environ.get(
                                         "INDUCTIVA_OUTPUT_DIR",
                                         "inductiva_output"))
_api_key = contextvars.ContextVar("INDUCTIVA_API_KEY",
                                  default=os.environ.get(
                                      "INDUCTIVA_API_KEY", None))

# Disable urllib3 warnings.
# TODO: Verify and fix the appearance of this warning.
urllib3_logger = logging.getLogger("urllib3.connectionpool")
urllib3_logger.setLevel(logging.CRITICAL)

__version__ = "0.8.4"


def set_output_dir(new_output_dir):
    """Sets the value of `inductiva._output_dir` to `new_output_dir`"""
    _output_dir.set(new_output_dir)


def get_output_dir():
    """Returns the value of inductiva._output_dir"""
    return _output_dir.get()


def _check_for_available_package_update():
    # pylint: disable=import-outside-toplevel
    from .localization import translator as __
    import urllib3
    import json

    new_version = __version__

    try:
        http = urllib3.PoolManager()
        resp = http.request("GET",
                            "https://pypi.org/pypi/inductiva/json",
                            headers={"Accept": "application/json"})
        json_response = json.loads(resp.data)
        new_version = json_response["info"]["version"]

    except Exception as ex:  # pylint: disable=broad-exception-caught
        logging.warning(__("failed-update-check", ex), exc_info=True)

    if new_version != __version__:
        msg = __("upgrade-available", new_version, __version__)
        print(msg, file=sys.stderr)


def _set_key_and_check_version():
    """Sets the api key and checks if it is valid."""
    if not utils.format_utils.getenv_bool("GITHUB_ACTIONS", False):
        set_api_key(get_api_key())

    # Perform version check only on first invocation
    if not hasattr(_set_key_and_check_version, "version_checked"):
        compare_client_and_backend_versions(__version__)
        _set_key_and_check_version.version_checked = True


_check_for_available_package_update()


def compare_client_and_backend_versions(client_version: str):
    """ Compares the provided client version 7with the backend API version.

    Sends a GET request to the backend API's version comparison endpoint
    with the client version as a parameter. Evaluates the response to
    determine if the client version is compatible with the backend version.
    Raises exceptions for communication issues or incompatibility.

    Parameters:
    - client_version (str): The version of the client to be compared with the
                            backend version.

    Raises:
    - RuntimeError: If the API cannot be reached, or if the client version is
      incompatible with the backend version, or for other general failures.
    """
    api_config = Configuration(host=api_url)

    with ApiClient(api_config) as client:
        api_instance = VersionApi(client)
        query_params = {"client_version": client_version}

        try:
            api_instance.compare_client_and_backend_versions(
                query_params=query_params)

        except (MaxRetryError, NewConnectionError) as exc:
            raise RuntimeError(
                "Failed to reach the API. "
                "Please check your connection and try again.") from exc

        except ApiException as e:
            if e.status == 406:
                raise RuntimeError(
                    f"Client version {client_version} is not compatible "
                    f"with API version {e.headers['version']}.\n"
                    "Please update the client version.") from e
            raise RuntimeError(e) from e

        except Exception as e:
            raise RuntimeError(
                f"Failed to compare client and API versions. {e}") from e


def set_api_key(api_key):
    """Sets the value of `inductiva._api_key` to `api_key"""
    if api_key is None:
        # pylint: disable=line-too-long
        raise ValueError(
            "No API Key specified. "
            "Please set the INDUCTIVA_API_KEY environment variable.\n"
            "More infomation at:"
            "https://console.inductiva.ai/")

    _api_key.set(api_key)


def get_api_key():
    """Returns the value of inductiva._api_key"""
    return _api_key.get()


def _supports_ansi():
    """Checks if we support ansi formatting for colors and bolds"""
    user_disable_ansi = utils.format_utils.getenv_bool("INDUCTIVA_DISABLE_ANSI",
                                                       False)
    if sys.platform.startswith("win"):
        return "TERM" in os.environ and os.environ[
            "TERM"] == "xterm" and not user_disable_ansi
    return hasattr(sys.stdout,
                   "isatty") and sys.stdout.isatty() and not user_disable_ansi


def _check_user_info():

    if utils.format_utils.getenv_bool("GITHUB_ACTIONS",False) or \
       utils.format_utils.getenv_bool("INDUCTIVA_DISABLE_IMPORT_INFO",False):
        return

    if not logs.is_cli():
        # Determine if we are importing from cli or script file
        # and only print info if called from a script file
        get_info(None, sys.stdout)


_ansi_enabled = _supports_ansi()

_set_key_and_check_version()

_check_user_info()
