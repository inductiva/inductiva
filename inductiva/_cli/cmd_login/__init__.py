"""Register CLI commands for Login in Inductiva API."""
import subprocess
import argparse
import platform

from .. import utils
from ... import constants


def api_login(unused_args):
    """Open the Inductiva API console to perform login."""

    url = constants.CONSOLE_URL

    system = platform.system()
    if system == "Darwin":
        cmd = ["open", url]
    elif system == "Linux":
        cmd = ["xdg-open", url]
    elif system == "Windows":
        cmd = ["cmd", "/c", "start", url.replace("&", "^&")]
    else:
        raise RuntimeError(f"Unsupported system: {system}")
    subprocess.run(cmd, check=True)


def register(root_parser):
    parser = root_parser.add_parser(
        "login",
        help="Login to get an Inductiva API Key.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.description = (
        "Login into the Inductiva console to register in the platform and "
        "obtain\nyour API key. Currently, this method requires an Google "
        "email.\nFurther, you can login to generate a new API key when needed."
        "\n\nTo use the API set the key as an environment variable:\n "
        "\texport INDUCTIVA_API_KEY=<your_api_key>")

    utils.show_help_msg(parser)

    parser.set_defaults(func=api_login)
