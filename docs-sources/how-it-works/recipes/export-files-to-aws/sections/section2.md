# Advanced Options

## Rename the Exported File in AWS S3
Specify a new filename during export:

```
inductiva storage export <file_path> --bucket-name <bucket_name> --file-name <new_file_name>
```

*Example:*

```
inductiva storage export simulation_results.tar.gz --bucket-name my-s3-bucket --file-name results_backup.tar.gz
```

## Optimize Uploads with Chunking Size
Adjust the minimum part size (in MB) for chunked uploads:

```
inductiva storage export <file_path> --bucket-name <bucket_name> --part-size <size>
```

*Example*:

```
inductiva storage export large_dataset.csv --bucket-name my-s3-bucket --part-size 100
```

This sets each upload chunk to be at least 100 MB, which can improve performance for large files.