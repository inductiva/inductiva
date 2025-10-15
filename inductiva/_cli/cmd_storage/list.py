"""List the storage contents via CLI."""

import argparse
import textwrap

from inductiva import storage


def listdir(args):
    """List the user's remote storage contents."""
    max_results = None if args.all else args.max_results

    storage.listdir(
        path=args.path,
        region=args.region,
        max_results=max_results,
        order_by=args.order_by,
        sort_order=args.sort_order,
    )


def register(parser):
    """Register the list user's remote storage command."""

    subparser = parser.add_parser(
        "list",
        aliases=["ls"],
        help="List remote storage contents.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    subparser.description = (
        "The `inductiva storage list` command provides an overview"
        " of your data on Inductiva remote storage.\n"
        "It lists all files and folders in a specified path, allowing you to "
        "control the maximum number of results,\n"
        "the ordering criteria, and the sorting order.\n")

    subparser.add_argument("path", default="/", type=str, nargs="?")

    subparser.add_argument("-m", "--max-results", default=10, type=int)

    subparser.add_argument(
        "-o",
        "--order-by",
        default="creation_time",
        type=str,
        choices=["creation_time", "size"],
        help="Order by creation_time or size.",
    )

    subparser.add_argument(
        "-s",
        "--sort-order",
        default="desc",
        type=str,
        choices=["desc", "asc"],
        help="Sorting order (desc or asc).",
    )

    subparser.add_argument(
        "--all",
        action="store_true",
        help="List all results, ignoring --max-results.",
    )

    subparser.add_argument(
        "-r",
        "--region",
        default=None,
        type=str,
        help=("Filter by region. Specify 'all' to list all regions. If not "
              "specified, the contents of the user's default region are "
              "returned."),
    )

    subparser.epilog = textwrap.dedent("""\
        examples:
            # List the 10 largest folders sorted by size
            $ inductiva storage list --max-results 10 --order-by size --sort-order desc

            NAME                             SIZE           CREATION TIME      PROVIDER     REGION
            0bet8jrpp2gz974n42nsd9n2p/       56.11 MB       06 Feb, 11:32:29   GCP          europe-west1
            05ujj5m0ytdkckxwk1tq1b5io/       27.93 MB       08 Feb, 09:19:44   GCP          europe-west1
            6a2h1wnxywpea8jfoxiikdjf7/       26.49 MB       07 Feb, 13:47:03   GCP          europe-west1
            f8joznwc9xf9a4nypcaei6v2s/       12.79 MB       07 Feb, 09:16:55   GCP          europe-west1
            dpq2cv6b5f9p1c77nc8anjo10/       12.00 MB       08 Feb, 09:39:31   GCP          europe-west1
            r4kerxf4b53krgn0s3fyece3b/       11.92 MB       07 Feb, 11:47:48   GCP          europe-west1
            j9qzrpiohgt7x97od3tw4wccd/       11.74 MB       07 Feb, 11:47:46   GCP          europe-west1
            iqi71gonoacfj7fknox3rvnq2/       11.52 MB       07 Feb, 11:47:45   GCP          europe-west1
            dxmnxdrfrv84pfbzbvm9v0dat/       11.43 MB       07 Feb, 11:47:43   GCP          europe-west1
            bgtwgnnyq5qa5hecegzdx6okr/       11.36 MB       07 Feb, 11:47:40   GCP          europe-west1

            Total storage size used:
                Volume: 5.31 GB
                Cost: 0.099 US$/month

            Listed 10 folder(s). Ordered by size.
            Use --max-results/-m to control the number of results displayed.

            You have storage in the following regions: europe-west1
    """)

    subparser.set_defaults(func=listdir)
