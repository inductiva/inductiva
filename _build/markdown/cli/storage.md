# storage

## inductiva storage - CLI interface

```default
inductiva storage [-h] {bucket,download,dl,export,list,ls,remove} ...
```

Remote storage management utilities. 

Use the `inductiva storage` command to manage your data in Inductiva’s remote storage. It supports a range of operations including listing files, downloading and deleting items, and calculating total storage usage and cost.

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva storage bucket

```default
inductiva storage bucket [-h] {list} ...
```

Storage bucket management utilities.

The `inductiva storage bucket` command provides tools for managing with cloud storage buckets. Use it to manage cloud storage buckets linked to your account.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva storage bucket list (ls)

```default
inductiva storage bucket list [-h]
```

The `inductiva storage bucket list` command retrieves and displays your remote storage buckets available through the Inductiva API. This includes information such as the bucket name, region, provider, and whether it is internally managed by Inductiva.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva storage download (dl)

```default
inductiva storage download [-h] [--decompress] [-r REGION] remote_path [local_dir]
```

The `inductiva storage download` command allows you to download a file or folder from your remote storage, maintaining the original storage path structure, to your local machine. Specify the remote path, the local destination (optional), and whether to decompress the content after downloading (default: enabled).

#### Positional Arguments

* [**`remote_path`**]() - The path to the file or folder in remote storage to download. (default: `None`)
* [**`local_dir`**]() - The local directory where the downloaded content will be saved. Defaults to the current working directory if not specified. (default: ``)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--decompress`**]() - Decompress the downloaded file or folder if it is compressed.
* [**`-r`**]() `REGION`, [**`--region`**]() `REGION` - Storage region of remote files. If not specified, the user’s default region is assumed.

#### Examples

```bash

# Download and decompress an archived file
$ inductiva storage download my_data/archive.zip --decompress
```

---

### inductiva storage export

```default
inductiva storage export [-h] [-r REGION] [--export-to {aws-s3}] [--file-name FILE_NAME]
                         --bucket-name BUCKET_NAME [--part-size PART_SIZE]
                         path_to_export
```

The `inductiva storage export` command lets you export data from your Inductiva remote storage to an external cloud provider, such as AWS S3.

To export to AWS S3, follow these steps:
1. Install `inductiva` with `pip install inductiva[aws]`.
2. Configure your AWS credentials using `aws configure`.
3. Ensure the target S3 bucket exists and you have write permissions.

#### Positional Arguments

* [**`path_to_export`**]() - File or folder path in Inductiva remote storage to export. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-r`**]() `REGION`, [**`--region`**]() `REGION` - Storage region of Inductiva remote files. If not specified, the user’s default region is assumed.
* [**`--export-to`**]() `EXPORT_TO` - External cloud service to export your data to. (default: `aws-s3`)
* [**`--file-name`**]() `FILE_NAME` - Name to assign to the file being saved. (default: `None`)
* [**`--bucket-name`**]() `BUCKET_NAME` - Name of the external cloud bucket where the file will be saved. (default: `None`)
* [**`--part-size`**]() `PART_SIZE` - Specify the size (in MB) of each part in the multipartupload. The default is 128 MB. For example, specify 50 for 50 MB.

#### Examples

```bash

# Export a file to AWS S3
inductiva storage export --export-to aws-s3 --bucket-name my-bucket my_data/file1.txt
```

---

### inductiva storage list (ls)

```default
inductiva storage list [-h] [-m MAX_RESULTS] [-o {creation_time,size}] [-s {desc,asc}]
                       [--all] [-r REGION]
                       [path]
```

The `inductiva storage list` command provides an overview of your data on Inductiva remote storage. It lists all files and folders in a specified path, allowing you to control the maximum number of results, the ordering criteria, and the sorting order.

#### Positional Arguments

* [**`path`**]() (default: `/`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-m`**]() `MAX_RESULTS`, [**`--max-results`**]() `MAX_RESULTS` (default: `10`)
* [**`-o`**]() `ORDER_BY`, [**`--order-by`**]() `ORDER_BY` - Order by creation_time or size. (default: `creation_time`)
* [**`-s`**]() `SORT_ORDER`, [**`--sort-order`**]() `SORT_ORDER` - Sorting order (desc or asc). (default: `desc`)
* [**`--all`**]() - List all results, ignoring –max-results.
* [**`-r`**]() `REGION`, [**`--region`**]() `REGION` - Filter by region. Specify `'all'` to list all regions. If not specified, the contents of the user’s default region are returned.

#### Examples

```bash

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
```

---

### inductiva storage remove (rm)

```default
inductiva storage remove [-h] [-y] [-a] [-r REGION] [paths ...]
```

The `inductiva storage remove` command deletes specific files orfolders from your Inductiva remote storage. Use it with caution. This action is irreversible and will permanently remove the selected files or directories from your remote storage.

#### Positional Arguments

* [**`paths`**]() - Remote path(s) to remove from storage. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-y`**](), [**`--yes`**]() - Sets any confirmation values to `"yes"` automatically. Users will not be asked for confirmation to remove path(s) from remote storage.
* [**`-a`**](), [**`--all`**]() - Remove all data from all remote storage regions.
* [**`-r`**]() `REGION`, [**`--region`**]() `REGION` - Storage region of Inductiva remote files. If not specified, the user’s default region is assumed. If –all flag is also specified this value will be ignored.

#### Examples

```bash

$ inductiva storage remove 0bet8jrpp2gz974n42nsd9n2p/
You are about to remove the following paths from your remote storage space:
- 0bet8jrpp2gz974n42nsd9n2p/
Are you sure you want to proceed (y/[N])? y
Removing '0bet8jrpp2gz974n42nsd9n2p/' from remote storage...
Successfully removed '0bet8jrpp2gz974n42nsd9n2p/' from remote storage.
```

---

### inductiva storage size

```default
inductiva storage size [-h]
```

The `inductiva storage size` command calculates the total size of your data on the platform. It returns the total size in GB of all items in your storage across all storage regions.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva storage upload (ul)

```default
inductiva storage upload [-h] [-r REGION] local_path remote_dir
```

The `inductiva storage upload` command allows you to upload files or folders to your remote storage. Specify the local path and the remote destination.

Example:
    inductiva storage upload local/path/file_or_directory remote_dir

#### Positional Arguments

* [**`local_path`**]() - The local path to the file or folder to upload. (default: `None`)
* [**`remote_dir`**]() - The remote directory where the uploaded content will be stored. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-r`**]() `REGION`, [**`--region`**]() `REGION` - Storage region of remote files. If not specified, the user’s default region is assumed.

---
::docsbannersmall
::
