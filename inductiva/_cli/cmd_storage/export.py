"""Export the user's remote storage to another cloud."""

import argparse
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

    subparser.description = (
        "The `export` command allows you to export your data to another cloud, "
        "such as AWS S3.\n  To export to AWS S3, you need to:\n"
        "1. Install Inductiva with `pip install inductiva[aws]`.\n"
        "2. Configure your AWS credentials using `aws configure`.\n"
        "3. Ensure the target S3 bucket exists and you have write permissions.")
    subparser.add_argument(
        "path_to_export",
        type=str,
        help="Specify the path of the file to export.",
    )
    subparser.add_argument(
        "--export-to",
        default=storage.ExportDestination.AWS_S3,
        type=storage.ExportDestination,
        choices=list(storage.ExportDestination),
        help="Specify the export destination: aws-s3.",
    )

    subparser.add_argument(
        "--file-name",
        type=str,
        required=False,
        help="Specify the name to assign to the file being saved.",
    )

    subparser.add_argument(
        "--bucket-name",
        type=str,
        required=True,
        help="Bucket name where to save the file.",
    )

    subparser.add_argument(
        "--part-size",
        type=int,
        required=False,
        default=128,
        help=("Specify the size (in MB) of each part in the multipartupload."
              "The default is 128 MB. For example, specify 50 for 50 MB."),
    )

    subparser.set_defaults(func=export)
