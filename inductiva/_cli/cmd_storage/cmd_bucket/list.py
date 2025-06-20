"""List the storage contents via CLI."""

import argparse

import inductiva.client
import inductiva
from inductiva.utils import format_utils


def list_buckets(_):
    """List the user's buckets accessible via the Inductiva API."""
    api = inductiva.client.StorageApi(inductiva.api.get_client())
    response_list_buckets = api.list_buckets()
    response_default_bucket = api.get_default_bucket()

    buckets = [schema.model_dump() for schema in response_list_buckets]
    default_bucket = response_default_bucket.model_dump()
    for bucket in buckets:
        bucket["is_default"] = "yes" if bucket == default_bucket else "no"
        bucket["is_internal"] = "yes" if bucket["is_internal"] else "no"

    cols = list(buckets[0].keys())
    rows = [list(bucket.values()) for bucket in buckets]

    ansi_formatter = format_utils.get_ansi_formatter()
    headers = [lambda x: ansi_formatter(x.upper(), format_utils.Emphasis.BOLD)]

    table = format_utils.get_tabular_str(rows, cols, header_formatters=headers)
    print(table)


def register(parser):
    """
    Register the `inductiva storage bucket list` command for listing
    user buckets.
    """

    subparser = parser.add_parser(
        "list",
        aliases=["ls"],
        help="List all your storage buckets accessible via the Inductiva API.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva storage bucket list` command retrieves and displays "
        "your remote storage buckets available through the Inductiva API.\n"
        "This includes information such as the bucket name, region, provider, "
        "and whether it is internally managed by Inductiva.\n"
        "Example:\n"
        "  inductiva storage bucket list\n")

    subparser.set_defaults(func=list_buckets)
