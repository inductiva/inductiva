"""Methods to interact with the user info on the API."""
import json
from typing import Any, Dict, List

import inductiva
from inductiva import api
from inductiva.client import ApiClient
from inductiva.client.apis.tags.users_api import UsersApi


def _fetch_quotas_from_api() -> List[Dict]:
    """Get information about a user's quotas.
    """
    api_config = api.validate_api_key(inductiva.api_key)

    with ApiClient(api_config) as client:
        api_instance = UsersApi(client)

        resp = api_instance.get_user_quotas().response
        quotas = json.loads(resp.data.decode("utf-8"))

        return quotas


def get_quotas() -> Dict[str, Dict[str, Any]]:
    """Get the user quotas.

    This function gets a dict with the user quotas.

    Returns:
        Dict with the user quotas.
    """
    quotas = _fetch_quotas_from_api()

    quotas_dict = {quota.pop("name"): quota for quota in quotas}

    return quotas_dict
