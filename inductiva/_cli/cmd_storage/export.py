"""Export the user's remote storage to another cloud."""

import argparse
import textwrap
from inductiva.storage import storage


def export(args):
    """Export the user's remote storage to another cloud."""
    storage.export(
        path_to_export=args.path_to_export,
        export_to=args.export_to,
        bucket_name=args.bucket_name,
        file_name=args.file_name,
        part_size=args.part_size,
    )


def register(parser):
    """Register the export user's storage command."""

    subparser = parser.add_parser(
        "export",
        help=("Copy files from Inductiva's storage to an external cloud "
              "bucket (e.g. AWS S3)."),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    subparser.description = textwrap.dedent("""\
        The `inductiva storage export` command lets you export data from
        your Inductiva remote storage to an external cloud provider, such as
        AWS S3.

        To export to AWS S3, follow these steps:
        1. Install `inductiva` with `pip install inductiva[aws]`.
        2. Configure your AWS credentials using `aws configure`.
        3. Ensure the target S3 bucket exists and you have write
        permissions.
    """)

    subparser.add_argument(
        "path_to_export",
        type=str,
        help="File or folder path in Inductiva remote storage to export.",
    )
    subparser.add_argument(
        "--export-to",
        default=storage.ExportDestination.AWS_S3,
        type=storage.ExportDestination,
        choices=list(storage.ExportDestination),
        help="External cloud service to export your data to.",
    )

    subparser.add_argument(
        "--file-name",
        type=str,
        required=False,
        help="Name to assign to the file being saved.",
    )

    subparser.add_argument(
        "--bucket-name",
        type=str,
        required=True,
        help="Name of the external cloud bucket where the file will be saved.",
    )

    subparser.add_argument(
        "--part-size",
        type=int,
        required=False,
        default=128,
        help=("Specify the size (in MB) of each part in the multipartupload. "
              "The default is 128 MB. For example, specify 50 for 50 MB."),
    )

    subparser.epilog = textwrap.dedent("""\
        examples:
            # Export a file to AWS S3
            inductiva storage export --export-to aws-s3 --bucket-name my-bucket my_data/file1.txt
    """)

    subparser.set_defaults(func=export)
