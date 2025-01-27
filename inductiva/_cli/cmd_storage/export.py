"""Export the user's remote storage to another cloud."""

import argparse
import pathlib
import inductiva
from enum import Enum
from inductiva.client.apis.tags import storage_api
import math
import logging

_boto3_imported = True
try:
    import boto3

    logging.getLogger("botocore").setLevel(logging.WARNING)
except ImportError:
    _boto3_imported = False


class ExportDestination(Enum):
    AWS_S3 = "aws-s3"

    def __str__(self):
        return self.value


def initiate_multipart_upload(filename, bucket_name):
    """
    Initiate a multipart upload on S3 and return the UploadId.
    """
    s3_client = boto3.client("s3")
    response = s3_client.create_multipart_upload(
        Bucket=bucket_name,
        Key=filename,
    )
    return response["UploadId"]


def generate_presigned_url(upload_id, part_number, filename, bucket_name):
    """
    Generate a presigned URL for uploading a part to S3.
    """
    s3_client = boto3.client("s3")
    method_parameters = {
        "Bucket": bucket_name,
        "Key": filename,
        "PartNumber": part_number,
        "UploadId": upload_id,
    }

    return s3_client.generate_presigned_url(
        "upload_part",
        Params=method_parameters,
        ExpiresIn=3600,
    )


def generate_complete_multipart_upload_signed_url(
    upload_id,
    filename,
    bucket_name,
):
    """
    Generate a presigned URL for completing the multipart upload.
    """
    s3_client = boto3.client("s3")

    signed_url = s3_client.generate_presigned_url(
        ClientMethod="complete_multipart_upload",
        Params={
            "Bucket": bucket_name,
            "Key": filename,
            "UploadId": upload_id
        },
        HttpMethod="POST",
    )

    return signed_url


def get_file_size(file_path):
    api = storage_api.StorageApi(inductiva.api.get_client())

    contents = api.list_storage_contents({
        "path": file_path,
        "max_results": 2,
    }).body
    if len(contents) > 1:
        raise ValueError(f"Multiple files found at {file_path}. "
                         "Please specify a single file.")

    return list(contents.values())[0]["size_bytes"]


def get_multipart_parts(size, min_part_size=50 * 1024 * 1024):
    """
    Calculate the size of each part and the total number of parts

    The goal is to divide the data into parts `min_part_size` each:
    1. No more than 10,000 parts are created (maximum parts allowed by S3).
    2. The part size might be increased to avoid exceeding the part limit.

    Args:
        size (int): The total size of the file to be uploaded, in bytes.
        min_part_size (int): The minimum size of each part, in bytes.

    Returns:
        tuple: (part_size, part_count)
            - part_size (int): The size of each part in bytes.
            - part_count (int): The total number of parts.
    """
    max_parts = 10000

    if size <= min_part_size:
        return size, 1

    # Calculate the part size based on the smaller of two values:
    # - At least `min_part_size`
    # - Maximum size to ensure no more than 10,000 parts (size // 10000)
    max_allowed_part_size = size // max_parts
    part_size = max(min_part_size, max_allowed_part_size)

    part_count = math.ceil(size / part_size)

    return part_size, part_count


def multipart_upload(
    path,
    parts_size,
    upload_parts,
    complete_multipart_url,
):
    """
    Perform the multipart upload using the server.
    """
    api = storage_api.StorageApi(inductiva.api.get_client())

    api.export_multipart_files(
        body={
            "path": path,
            "parts_size": parts_size,
            "upload_parts": upload_parts,
            "complete_multipart_url": complete_multipart_url,
        })


def export_to_aws_s3(path_to_export, min_part_size_mb, filename, bucket_name):
    if not _boto3_imported:
        print("boto3 is not installed. Please run "
              "'pip install inductiva[aws]' to install it.")
        return
    try:
        boto3.client("sts").get_caller_identity()
    except Exception:  # noqa: BLE001
        print("AWS credentials not found. "
              "Please set your AWS credentials with 'aws configure'")
        return

    # Step 1: Get the file size
    file_size = get_file_size(path_to_export)

    # Step 2: Calculate the part size and count
    parts_size, parts_count = get_multipart_parts(
        file_size,
        min_part_size=min_part_size_mb * 1024 * 1024,
    )

    # Step 3: Initiate the multipart upload on aws
    upload_id = initiate_multipart_upload(filename, bucket_name)

    # Step 4: Generate presigned URLs for each part
    upload_parts = []
    for part_number in range(1, parts_count + 1):
        presigned_url = generate_presigned_url(upload_id, part_number, filename,
                                               bucket_name)
        upload_parts.append({
            "part_number": part_number,
            "part_url": presigned_url
        })

    # Step 5: Generate the complete multipart upload signed URL
    complete_multipart_url = generate_complete_multipart_upload_signed_url(
        upload_id,
        filename,
        bucket_name,
    )

    # Step 6: Ask the server to perform the multipart upload
    multipart_upload(
        path_to_export,
        parts_size,
        upload_parts,
        complete_multipart_url,
    )
    print(
        "Export is being done by inductiva server. You can close the terminal.")


def export(args):
    if not args.file_name_to_save:
        filename = pathlib.Path(args.path_to_export).name
    else:
        filename = args.file_name_to_save

    if args.export_to == ExportDestination.AWS_S3:
        print(f"Exporting {args.path_to_export} to {args.bucket_name}...")
        export_to_aws_s3(
            args.path_to_export,
            args.min_part_size_MB,
            filename,
            args.bucket_name,
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
        default=ExportDestination.AWS_S3,
        type=ExportDestination,
        choices=list(ExportDestination),
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
