# FAQ

## 1. Do I need to create an AWS S3 bucket before exporting files?

Yes. Make sure your bucket already exists and that you have write permissions before exporting.

## 2. How can I verify that my file was successfully uploaded?

You can veridy using the AWS CLI:

```
aws s3 ls s3://<bucket_name>/
```

Or by visiting **S3 → Your Bucket → Objects** in the AWS Console.

## 3. Can I get a notification when the export finishes?

Yes — set up [S3 Event Notifications](https://docs.aws.amazon.com/AmazonS3/latest/userguide/EventNotifications.html) in AWS to receive alerts.

## 4. Is there a size limit for exported files?

Yes. Files up to **5 TB** are supported.

You can check the file size in the [Web Console](https://console.inductiva.ai/storage/storage) or in the **API** directly.

## 5. How can Inductiva export files if it doesn't have my AWS credentials?

The Inductiva client uses your local AWS credentials to [generate pre-signed URLs](https://docs.aws.amazon.com/AmazonS3/latest/userguide/ShareObjectPreSignedURL.html). These URLs are securely sent to Inductiva servers, allowing them to upload the file without ever accessing your credentials.

[AWS Multipart Upload](https://docs.aws.amazon.com/AmazonS3/latest/userguide/mpuoverview.html) is also used.

## 6. Can I use credentials limited to a specific bucket?

Yes, you can configure AWS credentials to restrict access to a specific S3 bucket. To do this, create an IAM policy that grants the `s3:PutObject` and `s3:ListBucket` (for uploading files) permissions specifically for that bucket.

## 7. How much does exporting to AWS S3 cost?

The data transfer costs associated with exporting files from Inductiva to AWS S3 will be included under your [storage costs](https://console.inductiva.ai/storage/storage).

## 8. Can I export multiple files at once?

Not currently. Each export command handles one file at a time.

## 9. Can I export files to cloud providers other than AWS?

Currently, only AWS S3 is supported.
