"""List the user costs via CLI."""
from decimal import Decimal
from typing import TextIO
import argparse
import sys

from inductiva import users, _cli
from inductiva.utils import format_utils


def get_costs(args, fout: TextIO = sys.stdout):
    """ Lists the user's costs.

    Lists all costs for the user in the specified time period.
    """
    start_year = args.start_year
    start_month = args.start_month
    end_year = args.end_year
    end_month = args.end_month
    user_costs = users.get_costs(start_year=start_year,
                                 start_month=start_month,
                                 end_year=end_year,
                                 end_month=end_month)

    for element in user_costs:
        compute_cost = Decimal(element["components"]["compute"])
        storage_cost = Decimal(element["components"]["storage"])
        total_cost = Decimal(element["total"])

        print("", file=fout)
        print("Month: ", element["month"], file=fout)
        print("\tType: ", element["type"], file=fout)
        print("\tTotal: ",
              format_utils.currency_formatter(total_cost),
              file=fout)
        print("\tComponents:", file=fout)
        print("\t\tCompute: ",
              format_utils.currency_formatter(compute_cost),
              file=fout)
        print("\t\tStorage: ",
              format_utils.currency_formatter(storage_cost),
              file=fout)

        warning = str(element.get("warning"))
        if "None" not in warning:
            print("\tWarning: ", element["warning"], file=fout)
    return 0


def register(parser):
    """Register the user costs command."""

    subparser = parser.add_parser("costs",
                                  help="List the user's costs in the specified "
                                  "time period.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva user costs` command provides "
        "an overview of your costs for a specified period.")

    subparser.add_argument(
        "--start-year",
        "-sy",
        type=int,
        required=True,
        help="The start year of the period to list the costs for.",
    )
    subparser.add_argument(
        "--start-month",
        "-sm",
        type=int,
        required=True,
        help="The start month of the period to list the costs for.",
    )
    subparser.add_argument(
        "--end-year",
        "-ey",
        type=int,
        required=True,
        help="The end year of the period to list the costs for.",
    )
    subparser.add_argument(
        "--end-month",
        "-em",
        type=int,
        required=True,
        help="The end month of the period to list the costs for.",
    )

    _cli.utils.add_watch_argument(subparser)

    subparser.set_defaults(func=get_costs)
