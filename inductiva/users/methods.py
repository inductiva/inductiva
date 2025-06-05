"""Methods to interact with the user info on the API."""
from typing import Dict, Optional, List

from inductiva import api

import inductiva.client
import inductiva.client.models


def _fetch_quotas_from_api() -> Dict[str, inductiva.client.models.Quota]:
    """Get information about a user's quotas.
    """

    with api.get_client() as client:
        api_instance = inductiva.client.UsersApi(client)

        return api_instance.get_user_quotas()


def get_quotas() -> Dict[str, inductiva.client.models.Quota]:
    """Get the user quotas.

    This function gets a dict with the user quotas.

    Returns:
        Dict with the user quotas.
    """
    return _fetch_quotas_from_api()


def get_info() -> inductiva.client.models.User:
    """Get the user information.

    This funtion gets the user information, including the user's name, email,
    username, plan, programs, and total available credits.

    Returns:
        Dict with the user information.
    """

    with api.get_client() as client:
        api_instance = inductiva.client.UsersApi(client)
        user = api_instance.get_user_info()
    return user


def get_costs(
    start_year: int,
    start_month: int,
    end_year: Optional[int] = None,
    end_month: Optional[int] = None,
) -> List[inductiva.client.models.CostDetail]:
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
        api_instance = inductiva.client.UsersApi(client)

        response = api_instance.get_user_costs(
            start_year=start_year,
            start_month=start_month,
            end_year=end_year,
            end_month=end_month,
        )

    return response.costs
