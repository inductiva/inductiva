# Export my files from Inductiva to AWS S3

Inductiva now makes it easy to export files from your [Personal Remote Storage](https://docs.inductiva.ai/en/latest/cli/access-storage.html) to AWS S3, streamlining your data management process.

With this feature, you can seamlessly copy your files to your AWS S3 bucket with just a few simple commands. Whether you’re backing up important results, sharing large datasets, or moving your files for further analysis, this integration with AWS S3 simplifies the export process.

## How to Use:

1. Install Inductiva with AWS Support

    Before exporting, ensure you have the necessary dependencies installed:

    `pip install inductiva[aws]`
    
    This step is required only once to enable AWS support.

2. Configure AWS Credentials

    To interact with AWS S3, set up your credentials using:

    `aws configure`

    You’ll be prompted to enter your AWS Access Key ID, Secret Access Key, Region, and Output format (e.g., JSON). Your AWS credentials are never sent to our server.

    For more details see [FAQ #5](#faqs).

3. Export Files to AWS S3 using Inductiva CLI

    Use the export command to transfer files from Inductiva’s storage to AWS S3 using `inductiva CLI`:

    `inductiva storage export <file_path> --bucket-name <bucket_name> --export-to aws-s3`

    Example

    `inductiva storage export <task_id>/output.zip --bucket-name my-s3-bucket --export-to aws-s3`

    This command uploads the output of the task to `my-s3-bucket` in AWS S3.

4. Export File to AWS S3 using Inductiva API

    You can also export files from Inductiva’s storage to AWS S3 programmatically using the Inductiva API. Below is a Python example demonstrating how to achieve this:

    ```python
    import inductiva

    inductiva.storage.export(
        path_to_export="<file_path>",
        export_to=inductiva.storage.ExportDestination.AWS_S3,
        bucket_name="my-s3-bucket",
    )
    ```

    This snippet uploads `<file_path>` which can be either a Task File or a [Remote File](https://tutorials.inductiva.ai/how_to/reuse-files.html) to `my-s3-bucket` in AWS S3.


## Advanced Options

- Rename the exported file in the AWS destination

    To save the file on AWS S3 with a different name while exporting:

    `inductiva storage export <file_path> --bucket-name <bucket_name> --file-name-to-save <new_file_name>`

    Example:

    `inductiva storage export simulation_results.tar.gz --bucket-name my-s3-bucket --file-name-to-save results_backup.tar.gz`

- Optimize Upload with Chunking size

    The export commands uploads files in chunks. You can specify the minimum part size (in MB) to potentially upload multiple parts in parallel:

    `inductiva storage export <file_path> --bucket-name <bucket_name> --min-part-size-MB <size>`

    Example:

    `inductiva storage export large_dataset.csv --bucket-name my-s3-bucket --min-part-size-MB 100`

    This configures each chunk to be at least 100MB, optimizing for large file transfers.


## FAQs


1. Do I need an existing AWS S3 bucket?

    Yes. Ensure your S3 bucket exists and that you have write permissions before exporting.


2. How can I verify that my file was successfully uploaded?

    You can check your AWS S3 bucket using:

    `aws s3 ls s3://<bucket_name>/`

    Or through the AWS Console under `S3 → Your Bucket → Objects.`

3. Can I get a notification when the export finishes?

    Yes, you can use [S3 Event Notications](https://docs.aws.amazon.com/AmazonS3/latest/userguide/EventNotifications.html) to trigger an action when a file is uploaded to the bucket.

4. Is there a size limit for exporting files?

    Currently the export limit is 5TB per file.

    You can check the file size in the [Web Console](https://console-dev.inductiva.ai/storage) or in the [API directly](https://docs.inductiva.ai/en/latest/cli/access-storage.html#list-storage-contents).

5. How can Inductiva export files to AWS S3 if it does not have my AWS credentials?

    The inductiva client generates [aws pre-signed urls](https://docs.aws.amazon.com/AmazonS3/latest/userguide/ShareObjectPreSignedURL.html) with the AWS credentials which are sent to the server. The server then uploads the files to those pre-signed urls without requiring the AWS credentials. [AWS Multipart Upload](https://docs.aws.amazon.com/AmazonS3/latest/userguide/mpuoverview.html) is also used.

    This means that the Inductiva server never receives your AWS credentials.

6. Can I use AWS credentials with access limited to a specific bucket?

    Yes, you can configure the AWS credentials to provide restricted access to a specific S3 bucket. To do this, create an IAM policy that grants the `s3:PutObject` permission (for uploading files) specifically for that bucket. 

7. How much does exporting to AWS S3 cost?

    The data transfer costs associated with exporting files from Inductiva to AWS S3 will be included under your [storage costs](https://console-dev.inductiva.ai/account/costs).

7. Can I export multiple files at once?

    Currently, you need to run the export command for each file.

8. Can I export files to cloud providers other than AWS?

    At the moment, only AWS S3 is supported.
