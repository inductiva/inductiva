# storage

The `inductiva storage` command allows you to manage your remote storage on the Inductiva platform.
You can list, download, remove, and calculate the total storage usage of your stored data.

## Usage

```sh
inductiva storage [-h] {download,dl,export,list,ls,remove,rm,size} ...
```

### Description
Inductiva provides a remote storage system where users can store and manage
files related to their tasks. This command lets you interact with your stored data efficiently.

## Options

- **`-h, --help`** → Show help message and exit.

## Available Subcommands

### `list (ls)`
List all files and folders in your remote storage.

```sh
inductiva storage list
```

or using the alias:

```sh
inductiva storage ls
```

### `download (dl)`
Download a file or folder from remote storage to your local machine.

```sh
inductiva storage download [-h] [--decompress] remote_path [local_dir]
```

#### Options for `download`:
- **`-h, --help`** → Show help message and exit.
- **`--decompress`** → Automatically decompress the downloaded file if it is archived.
- **`remote_path`** → Path to the file or folder in remote storage (required).
- **`local_dir`** → (Optional) Destination directory for the downloaded file or folder.

or using the alias:

```sh
inductiva storage dl <file_or_folder_path>
```

If `remote_path` is not provided, the command will return an error:

```sh
inductiva storage download: error: the following arguments are required: remote_path
```

### `export`
Copy files from Inductiva's storage to an external cloud bucket (e.g., AWS S3).

```sh
inductiva storage export [-h] [--export-to {aws-s3}] [--file-name FILE_NAME] --bucket-name BUCKET_NAME [--part-size PART_SIZE] path_to_export
```

#### Options for `export`:
- **`-h, --help`** → Show help message and exit.
- **`--export-to {aws-s3}`** → Specify the external cloud service (currently supports AWS S3).
- **`--file-name FILE_NAME`** → (Optional) Rename the file before exporting.
- **`--bucket-name BUCKET_NAME`** → Name of the external cloud bucket (required).
- **`--part-size PART_SIZE`** → (Optional) Define the part size for chunked uploads.
- **`path_to_export`** → Path of the file or folder in Inductiva storage that needs to be exported.

### `remove (rm)`
Delete specific files or folders from remote storage.

```sh
inductiva storage remove <file_or_folder_path>
```

or using the alias:

```sh
inductiva storage rm <file_or_folder_path>
```

Use this command with caution as it **permanently** 
deletes data from your remote storage. 

### `size`
Retrieve the total storage usage.

```sh
inductiva storage size
```

## Example Usage

### List all stored files:
```sh
inductiva storage list
```
List the 10 largest folders sorted by size.
```bash
$ inductiva storage list --max-results 10 --order-by size --sort-order desc

       NAME                             SIZE           CREATION TIME
       0bet8jrpp2gz974n42nsd9n2p/       56.11 MB       06 Feb, 11:32:29
       05ujj5m0ytdkckxwk1tq1b5io/       27.93 MB       08 Feb, 09:19:44
       6a2h1wnxywpea8jfoxiikdjf7/       26.49 MB       07 Feb, 13:47:03
       f8joznwc9xf9a4nypcaei6v2s/       12.79 MB       07 Feb, 09:16:55
       dpq2cv6b5f9p1c77nc8anjo10/       12.00 MB       08 Feb, 09:39:31
       r4kerxf4b53krgn0s3fyece3b/       11.92 MB       07 Feb, 11:47:48
       j9qzrpiohgt7x97od3tw4wccd/       11.74 MB       07 Feb, 11:47:46
       iqi71gonoacfj7fknox3rvnq2/       11.52 MB       07 Feb, 11:47:45
       dxmnxdrfrv84pfbzbvm9v0dat/       11.43 MB       07 Feb, 11:47:43
       bgtwgnnyq5qa5hecegzdx6okr/       11.36 MB       07 Feb, 11:47:40



### Download a file:
```sh
inductiva storage download my_data/file1.txt
```

### Download and decompress an archived file:
```sh
inductiva storage download my_data/archive.zip --decompress
```

### Remove a folder:
```sh
inductiva storage remove my_data/folder1
```

Remove the largest folder as listed above:

```sh
$ inductiva storage remove 0bet8jrpp2gz974n42nsd9n2p/
You are about to remove the following paths from your remote storage space:
  - 0bet8jrpp2gz974n42nsd9n2p/
Are you sure you want to proceed (y/[N])? y
Removing '0bet8jrpp2gz974n42nsd9n2p/' from remote storage...
Successfully removed '0bet8jrpp2gz974n42nsd9n2p/' from remote storage.
```

### Export a file to AWS S3:
```sh
inductiva storage export --export-to aws-s3 --bucket-name my-bucket my_data/file1.txt
```

### Check total storage usage:
```sh
inductiva storage size
```

## Need Help?
Run the following command for more details:

```sh
inductiva storage --help
```

