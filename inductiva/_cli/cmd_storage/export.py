"""Export the user's remote storage to another cloud."""

import argparse
from enum import Enum
from inductiva.storage import storage


def export(args):
    """Export the user's remote storage to another cloud."""
    storage.export(
        path_to_export=args.path_to_export,
        export_to=args.export_to,
        bucket_name=args.bucket_name,
        file_name_to_save=args.file_name_to_save,
        min_part_size_MB=args.min_part_size_MB,
    )


def register(parser):
    """Register the export user's storage command."""

    subparser = parser.add_parser(
        "export",
        help="Export the user's remote storage to another cloud.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    subparser.description = (
        "The `export` command allows you to export your data to another cloud, "
        "such as AWS S3.")
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
        "--file-name-to-save",
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
        "--min-part-size-MB",
        type=int,
        required=False,
        default=50,
        help=(
            "Specify the minimum size (in MB) of each part in the multipart "
            "upload. The default is 50 MB. For example, specify 50 for 50 MB."),
    )

    subparser.set_defaults(func=export)
