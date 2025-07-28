# inductiva storage [\[subcommands\]](#subcommands) [\[flags\]](#flags)

The `inductiva storage` command allows you to manage your remote storage on the Inductiva platform.

Inductiva provides a [remote storage](../../how-it-works/cloud-storage/index.md) system where users can store and manage
files related to their tasks. This command lets you list, download, remove, and calculate the total storage usage of your stored data.

````{eval-rst}
.. seealso::
   For complete API documentation, see the `Storage <https://inductiva.ai/guides/api-functions/api/inductiva.storage>`_ class documentation
````

## Subcommands

### `list (ls)` [\[flags\]](#flags-for-list)
List all files and folders in your remote storage.

This subcommand displays the contents of your remote storage at a given `path`. If **no path is provided**, it lists the contents of the **root directory (`/`)**. You can customize the number of results, the sort field, and the sorting order.

```sh
inductiva storage list
```

or using the alias:

```sh
inductiva storage ls
```

You can also specify a path:

```sh
inductiva storage ls [<PATH>]
```

List the 10 largest folders sorted by size:
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
```

<h4 id="flags-for-list">Flags</h4>

**`--max-results, -m`** (default: 10)

Limits the number of results returned. Ignored if `--all` is set.

---

**`--order-by, -o`** (default: `creation_time`)

Criteria to order results by. Options:

- `creation_time`: Sort by time the file or folder was created.
- `size`: Sort by file size (in bytes).

---

**`--sort-order, -s`** (default: `desc`)

Order direction of the listing. Options:

- `desc`: Descending.
- `asc`: Ascending.

---

**`--all`**

List all results.

### `download (dl)` [\[flags\]](#flags-for-download)
Download a file or folder from remote storage to your local machine.

By default, the downloaded content is saved to your **current working directory**. You can optionally specify a target local directory.

```sh
inductiva storage download <REMOTE_PATH> [<LOCAL_DIR>]
```
- `<REMOTE_PATH>` is the path to the file or folder in remote storage (required).
- `<LOCAL_DIR>` is the **optional** destination path on your local machine.

Download and decompress an archived file:

```sh
$ inductiva storage download my_data/archive.zip --decompress
```

<h4 id="flags-for-download">Flags</h4>

**`--decompress`**

Automatically decompress the downloaded file if it is archived.

### `export` [\[flags\]](#flags-for-export)
Copy files from Inductiva's storage to an external cloud bucket (e.g., AWS S3).

```sh
inductiva storage export <PATH_TO_EXPORT>
```

Export a file to AWS S3:

```sh
$ inductiva storage export --export-to aws-s3 --bucket-name my-bucket my_data/file1.txt
```

<h4 id="flags-for-download">Flags</h4>

**`--export-to {aws-s3}`**

Specify the external cloud service (currently supports AWS S3).

---

**`--file-name=<filename>`** (Optinal)

Rename the file before exporting.

---

**`--bucket-name=<bucket>`** (Required)

Name of the external cloud bucket.

---

**`--part-size=<size>`** (Optional)

Specify the size (in MB) of each part in the chunked upload. The default is 128 MB. For example, specify 50 for 50 MB.

### `remove (rm)` [\[flags\]](#flags-for-remove)
Delete specific files or folders from remote storage.

```sh
inductiva storage remove <FILE_OR_FOLDER_PATH>
```

or using the alias:

```sh
inductiva storage rm <FILE_OR_FOLDER_PATH>
```

Remove a specific folder:

```sh
$ inductiva storage remove 0bet8jrpp2gz974n42nsd9n2p/
  You are about to remove the following paths from your remote storage space:
    - 0bet8jrpp2gz974n42nsd9n2p/
  Are you sure you want to proceed (y/[N])? y
  Removing '0bet8jrpp2gz974n42nsd9n2p/' from remote storage...
  Successfully removed '0bet8jrpp2gz974n42nsd9n2p/' from remote storage.
```

Use this command with caution as it **permanently** 
deletes data from your remote storage. 

<h4 id="flags-for-remove">Flags</h4>

**`--yes, -y`**

Automatically answers **yes** to all confirmation prompts, allowing operations to proceed without interactive confirmation. Use with caution for irreversible actions.

---

**`--all, -a`**

Remove all data from remote storage.

### `size`
Retrieve the total storage usage.

```sh
inductiva storage size
```

## Flags
### `-h, --help`

Show help message and exit.

## Need Help?
Run the following command for more details:

```sh
inductiva storage --help
```


```{banner_small}
:origin: cli-storage
```
