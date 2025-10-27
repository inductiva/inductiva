"""Methods to interact with the user info on the API."""
from typing import Dict, Optional, List
from inductiva.utils import format_utils

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


def get_fees() -> inductiva.client.models.UserFees:
    """Get information about the user's fees."""
    with api.get_client() as client:
        api_instance = inductiva.client.UsersApi(client)
        fees = api_instance.get_user_fees()
    return fees


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


def get_task_orchestration_fee_warning():
    fees = get_fees()
    task_orchestration_fee = format_utils.currency_formatter(
        fees.task_orchestration_fee.amount)

    active_since_msg = ""
    if fees.task_orchestration_fee.active_since:
        task_orchestration_active_since = (
            format_utils.datetime_formatter_month_text(
                fees.task_orchestration_fee.active_since))

        active_since_msg = f" run from {task_orchestration_active_since}"

    return (
        f"Note: A per-run orchestration fee ({task_orchestration_fee}) "
        f"applies to tasks{active_since_msg}, in addition to the computation "
        "costs.\nLearn more about costs at: https://inductiva.ai/guides"
        "/how-it-works/basics/how-much-does-it-cost")
