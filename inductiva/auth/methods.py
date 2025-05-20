"""Methods to interact with the user info on the API."""
import argparse

from inductiva._cli.cmd_auth.login import login as login_cmd
from inductiva._cli.cmd_auth.logout import logout as logout_cmd


def login(private: bool = False) -> None:
    """Logs the user in to the Inductiva platform.

    This function handles user login by prompting the user to log in via the
      command line.

    Args:
        private:  If True, no api_key will be printed to the console.
    """

    args = argparse.Namespace(private=private)
    login_cmd(args)


def logout():
    """Logs the user out by removing the stored API key."""
    logout_cmd(None)
