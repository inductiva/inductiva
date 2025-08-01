"""List the user information via CLI."""
from collections import defaultdict
from datetime import datetime
import textwrap
from typing import TextIO, Dict
import argparse
import sys

from inductiva import users, _cli
import inductiva.client.models
from inductiva.client.exceptions import ApiException
from inductiva.users.methods import get_costs
from inductiva.utils import format_utils


def _quotas_to_dict(quotas: Dict[str, inductiva.client.models.Quota]):
    """Converts a list of quotas to two dictionaries.
    One for global quotas and one for instance quotas.
    """
    table_global = defaultdict(list)
    table_instance = defaultdict(list)
    for _, quota in quotas.items():
        current = quota.in_use
        max_allowed = quota.max_allowed
        unit = quota.unit
        current_str = (f"{current} {unit}" if current is not None else "N/A")
        max_allowed_str = (f"{max_allowed} {unit}"
                           if max_allowed is not None else "N/A")

        if quota.scope == "instance":
            table_instance[""].append(quota.label)
            table_instance["max allowed"].append(max_allowed_str)
        else:
            table_global[""].append(quota.label)
            table_global["current usage"].append(current_str)
            table_global["max allowed"].append(max_allowed_str)
    return table_global, table_instance


def _print_quotas(_, fout: TextIO = sys.stdout):
    """ Lists the user's quotas.

    Lists all the user's quotas and the quotas left for the user to use.
    """
    table, table_instance = _quotas_to_dict(users.get_quotas())

    emph_formatter = format_utils.get_ansi_formatter()
    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    table = format_utils.get_tabular_str(table,
                                         header_formatters=header_formatters)

    table_instance = format_utils.get_tabular_str(
        table_instance, header_formatters=header_formatters)
    print("■ Global User quotas", end="")
    print(table, file=fout)
    print("■ Instance User quotas", end="")
    print(table_instance, file=fout, end="")
    return 0


def _print_credits_summary(user: inductiva.client.models.User,
                           fout: TextIO = sys.stdout):
    """Prints the user's credits information."""
    total_available_credits = format_utils.currency_formatter(
        user.total_available_credits)
    print(f"■ Credits: {total_available_credits}", file=fout)


def _print_estimated_costs(fout: TextIO = sys.stdout):
    """Prints the user's estimated costs for the current month."""
    current_month = datetime.now().month
    current_year = datetime.now().year
    # This endpoint is not available for external users for now.
    # This try will be removed in the future.
    try:
        costs = get_costs(start_year=current_year, start_month=current_month)
        total_estimated_costs = format_utils.currency_formatter(costs[0].total)
        computation_estimated_costs = format_utils.currency_formatter(
            costs[0].components.compute)
        storage_estimated_costs = format_utils.currency_formatter(
            costs[0].components.storage)
        data_transfer_estimated_costs = format_utils.currency_formatter(
            costs[0].components.data_transfer)

        label_width = 12
        cost_width = 10

        print("■ Estimated Costs (current month):", file=fout)
        print(
            f"\t{'Computation:':<{label_width}} "
            f"{computation_estimated_costs:<{cost_width}}",
            file=fout)
        print(
            f"\t{'Storage:':<{label_width}} "
            f"{storage_estimated_costs:<{cost_width}}",
            file=fout)
        print(
            f"\t{'Data transfer:':<{label_width}} "
            f"{data_transfer_estimated_costs:<{cost_width}}",
            file=fout)
        print(
            f"\t{'Total:':<{label_width}} "
            f"{total_estimated_costs:<{cost_width}}",
            file=fout)
    except ApiException as _:
        return


def get_info(_, fout: TextIO = sys.stdout):
    """ Lists the user's credits.

    Lists all the user's credits and the credits left for the user to use.
    """
    user = users.get_info()

    organization = user.organization
    tier = user.tier
    username = user.username
    email = user.email
    name = user.name

    if organization:
        print(f"\nOrganization: {organization}\n", file=fout)

    print(f"Name: {name}", file=fout)
    print(f"Email: {email}", file=fout)
    print(f"Username: {username}", file=fout)

    print("", file=fout)
    print(f"■ Plan: {tier}", file=fout)
    print("", file=fout)

    _print_credits_summary(user, fout=fout)

    _print_estimated_costs(fout=fout)

    print("", file=fout)
    _print_quotas(_, fout)
    print("", file=fout)

    return 0


def register(parser):
    """Register the user credits command."""

    subparser = parser.add_parser("info",
                                  help="List the user's information.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = textwrap.dedent("""\
        The `inductiva user info` command displays details about your Inductiva
        account, including your user ID, email address, available credits,
        usage quotas, and other useful information.
    """)

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_info)
