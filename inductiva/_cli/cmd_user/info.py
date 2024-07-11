"""List the user information via CLI."""
from collections import defaultdict
from typing import TextIO
import argparse
import sys

from inductiva import users, _cli
from inductiva.utils import format_utils


def _quotas_to_dict(quotas):
    """Converts a list of quotas to two dictionaries.
    One for global quotas and one for instance quotas.
    """
    table_global = defaultdict(list)
    table_instance = defaultdict(list)
    for _, quota in quotas.items():
        current = quota["in_use"]
        max_allowed = quota["max_allowed"]
        unit = quota["unit"]
        current_str = (f"{current} {unit}" if current is not None else "N/A")
        max_allowed_str = (f"{max_allowed} {unit}"
                           if max_allowed is not None else "N/A")

        if quota["scope"] == "instance":
            table_instance[""].append(quota["label"])
            table_instance["max allowed"].append(max_allowed_str)
        else:
            table_global[""].append(quota["label"])
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


def _print_credits_summary(user_info, fout: TextIO = sys.stdout):
    """Prints the user's credits information."""
    total_available_credits = user_info["total_available_credits"]
    available_credits = user_info["tier"]["available_credits"]
    currency = user_info["credits_currency"]
    tier = user_info["tier"]["name"]

    print("■ Credits", file=fout)
    print("", file=fout)
    # pylint: disable=consider-using-f-string
    print("  {0:<25s} {1:10.2f} {2}".format(tier + " (tier)", available_credits,
                                            currency),
          file=fout)
    for program in user_info["campaigns"]:
        print("  {0:<25s} {1:10.2f} {2}".format(program["name"] + " (campaign)",
                                                program["available_credits"],
                                                currency),
              file=fout)

    print("  ----------------------------------------", file=fout)
    print("  {0:<25s} {1:10.2f} {2}".format("Total", total_available_credits,
                                            currency),
          file=fout)


def _campaigns_to_dict(campaigns):
    """Converts a list of campaigns to a dictionary."""
    table = defaultdict(list)
    for program in campaigns:
        table["name"].append(program["name"])
        table["enrollment date"].append(program["enrollment_date"])
        table["expiry date"].append(program["expiry_date"])
        table["available credits"].append(program["available_credits"])
        table["initial credits"].append(program["initial_credits"])
    return table


def get_info(_, fout: TextIO = sys.stdout):
    """ Lists the user's credits.

    Lists all the user's credits and the credits left for the user to use.
    """
    user_info = users.get_info()

    tier = user_info["tier"]["name"]
    username = user_info["username"]
    email = user_info["email"]
    name = user_info["name"] or ""

    print(f"Name: {name}", file=fout)
    print(f"Email: {email}", file=fout)
    print(f"Username: {username}", file=fout)

    print("", file=fout)
    print(f"■ Tier: {tier}", file=fout)
    print("", file=fout)

    _print_credits_summary(user_info, fout=fout)

    table = _campaigns_to_dict(user_info["campaigns"])

    emph_formatter = format_utils.get_ansi_formatter()

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    formatters = {
        "enrollment date": [format_utils.datetime_formatter_ymd_hm,],
        "expiry date": [format_utils.datetime_formatter_ymd_hm,]
    }

    table = format_utils.get_tabular_str(table,
                                         formatters=formatters,
                                         header_formatters=header_formatters)
    print("", file=fout)
    print("■ Active Campaigns", file=fout)
    if not user_info["campaigns"]:
        print("\n Not currently enrolled in any campaign.\n", file=fout)
    else:
        print(table, file=fout)

    _print_quotas(_, fout)

    print("", file=fout)
    return 0


def register(parser):
    """Register the user credits command."""

    subparser = parser.add_parser("info",
                                  help="List the user's information.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva user info` command provides "
        "an overview of your tier, campaigns and credits.\n"
        "It lists all your campaigns as well as the"
        "credits left for you to use.\n")

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_info)
