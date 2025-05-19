"""Methods to interact with the user info on the API."""
import argparse
import json
from typing import Any, Dict

from inductiva import api
from inductiva._cli.cmd_auth.login import login as login_cmd
from inductiva._cli.cmd_auth.logout import logout as logout_cmd

from inductiva.client.apis.tags.users_api import UsersApi


def _fetch_quotas_from_api() -> Dict[str, Dict[str, Any]]:
    """Get information about a user's quotas.
    """

    with api.get_client() as client:
        api_instance = UsersApi(client)

        resp = api_instance.get_user_quotas(skip_deserialization=True).response
        quotas = json.loads(resp.data.decode("utf-8"))

        return quotas


def get_quotas() -> Dict[str, Dict[str, Any]]:
    """Get the user quotas.

    This function gets a dict with the user quotas.

    Returns:
        Dict with the user quotas.
    """
    quotas = _fetch_quotas_from_api()

    return quotas


def get_info() -> Dict[str, Any]:
    """Get the user information.

    This funtion gets the user information, including the user's name, email,
    username, plan, programs, and total available credits.

    Returns:
        Dict with the user information.
    """

    with api.get_client() as client:
        api_instance = UsersApi(client)
        request = api_instance.get_user_info()
    return request.body


def get_costs(start_year: int,
              start_month: int,
              end_year: int = None,
              end_month: int = None) -> Dict[str, Any]:
    """Get the user costs.

    This function gets a dict with the user costs.

    Returns:
        Dict with the user costs.
    """
    if end_year and not end_month:
        raise ValueError("If end_year is provided, "
                         "end_month must also be provided.")
    if end_month and not end_year:
        raise ValueError("If end_month is provided, "
                         "end_year must also be provided.")

    with api.get_client() as client:
        api_instance = UsersApi(client)

        query_params = {"start_year": start_year, "start_month": start_month}

        #add end year and month if provided
        if end_year and end_month:
            query_params["end_year"] = end_year
            query_params["end_month"] = end_month

        request = api_instance.get_user_costs(query_params=query_params)
    return request.body["costs"]


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
