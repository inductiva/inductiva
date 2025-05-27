# Get Started in 3 Easy Steps

## 1. Install Inductiva with AWS Support

Before you begin, ensure you have the required dependencies installed:

```
pip install inductiva[aws]
```

> You only need to do this once to enable AWS support.

## 2. Configure Your AWS Credentials

Use the AWS CLI to set up your credentials for S3 access:

```
aws configure
```

You'll be prompted to enter:

- AWS Access Key ID
- AWS Secret Access Key
- Default region
- Output format (e.g., JSON)

> Your credentials are **never sent to Inductiva**. For more details, see [FAQ #5](https://inductiva.ai/guides/how-it-works/recipes/export-files-to-aws/sections/section3#how-can-inductiva-export-files-if-it-doesn-t-have-my-aws-credentials) below.

## 3. Export Files to AWS S3

### Option A: Using the Inductiva CLI

Use the `inductiva CLI` to export files from Inductiva storage to an AWS S3 bucket:

```
inductiva storage export <file_path> --bucket-name <bucket_name> --export-to aws-s3
```

*Example*:

```
inductiva storage export <task_id>/output.zip --bucket-name my-s3-bucket --export-to aws-s3
```

This uploads the task output to your specified AWS S3 bucket.

### Option B: Using the Inductiva API

Export programmatically with the Python API:

```python
import inductiva

inductiva.storage.export(
    path_to_export="<file_path>",
    export_to=inductiva.storage.ExportDestination.AWS_S3,
    bucket_name="my-s3-bucket",
)
```

This uploads `<file_path>`, which can be a task output or a remote file, to your AWS S3 bucket.
