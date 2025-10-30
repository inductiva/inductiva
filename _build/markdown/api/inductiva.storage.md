# storage

Methods to interact with the user storage resources.

### *class* ExportDestination(value)

Bases: `Enum`

#### AWS_S3 *= 'aws-s3'*

### *class* StorageOperation(api: StorageApi, id_)

Bases: `object`

Represents a storage operation running remotely via Inductiva API.

#### \_\_init_\_(api: StorageApi, id_)

#### *classmethod* from_api_response(api, response)

#### wait(poll_s: int = 2)

Wait for the operation to complete.

* **Parameters:**
  **poll_s** – Time in seconds between calls to the API to update
  the status of the operation.

### *class* ZipArchiveInfo(size: int, files: List[[ZipFileInfo](#inductiva.storage.storage.ZipFileInfo)])

Bases: `object`

Represents the total ZIP size and file contents of a ZIP archive.

#### \_\_init_\_(size: int, files: List[[ZipFileInfo](#inductiva.storage.storage.ZipFileInfo)]) → None

#### files *: List[[ZipFileInfo](#inductiva.storage.storage.ZipFileInfo)]*

#### size *: int*

### *class* ZipFileInfo(name: str, size: int | None, compressed_size: int | None, range_start: int | None, creation_time: datetime | None, compress_type: int | None)

Bases: `object`

Represents information about a file within a ZIP archive.

#### \_\_init_\_(name: str, size: int | None, compressed_size: int | None, range_start: int | None, creation_time: datetime | None, compress_type: int | None) → None

#### compress_type *: int | None*

#### compressed_size *: int | None*

#### creation_time *: datetime | None*

#### name *: str*

#### range_start *: int | None*

#### size *: int | None*

### *class* ZipFileRange(range_start: int, range_end: int)

Bases: `object`

Represents the byte range of a file within a ZIP archive.

#### \_\_init_\_(range_start: int, range_end: int) → None

#### range_end *: int*

#### range_start *: int*

### copy(source: str, target: str, region: str | None = None)

Copies a file or folder from a source path in storage to a target path.

* **Parameters:**
  * **source** (*str*) – The source path of the file or directory to copy.
  * **target** (*str*) – The destination path where the file or directory
    should be copied to.
  * **region** (*str* *,* *optional*) – The region of the remote storage. If not
    specified, the user’s default region is assumed.

### download(remote_path: str, local_dir: str = '', decompress: bool = True, region: str | None = None)

Downloads a file or folder from storage to a local directory, optionally
decompressing the contents.

* **Parameters:**
  * **remote_path** (*str*) – The path of the file or folder on the remote server
    to download.
  * **local_dir** (*str* *,* *optional*) – The local directory where the file or folder
    will be saved. Defaults to the current working directory.
  * **decompress** (*bool* *,* *optional*) – Whether to decompress the downloaded file
    or folder if it is compressed. Defaults to True.
  * **region** (*str* *,* *optional*) – The region of the remote storage. If not
    specified, the user’s default region is assumed.

### Examples

Download a folder from a remote server to the current directory:

```python
inductiva.storage.download(remote_path="/path/to/remote/folder/")
```

Download a file and save it to a local directory without decompressing:

```python
inductiva.storage.download(remote_path="/path/to/remote/file.zip",
                           local_dir="/local/directory",
                           decompress=False)
```

Download a file inside a zip archive:

```python
inductiva.storage.download(
    remote_path="/some_task_id/output.zip/stdout.txt"
)
```

#### NOTE
It is not possible to download folders that are inside zip archives.

### export(path_to_export: str, export_to: [ExportDestination](#inductiva.storage.storage.ExportDestination), bucket_name: str, file_name: str | None = None, part_size: int = 128, region: str | None = None)

### export_to_aws_s3(path_to_export, part_size, filename, bucket_name, region: str | None = None)

### get_file_range(path: str, zip_relative_path: str = '', filename: str = '', region: str | None = None) → [ZipFileRange](#inductiva.storage.storage.ZipFileRange)

Retrieve the byte range (start and end) of the compressed data for a
specific file inside a ZIP archive.

* **Parameters:**
  * **path** (*str*) – The full path to the ZIP archive.
  * **zip_relative_path** (*str* *,* *optional*) – The relative path inside the ZIP.
  * **filename** (*str*) – The name of the file inside the ZIP to get the range.
  * **region** (*str* *,* *optional*) – The region of the remote storage. If not
    specified, the user’s default region is assumed.
* **Returns:**
  The start and end byte offsets of the file
* **Return type:**
  [ZipFileRange](#inductiva.storage.storage.ZipFileRange)

### get_signed_urls(paths: List[str], operation: Literal['upload', 'download'], region: str | None = None) → List[str]

### get_space_used()

Returns the occupied storage size in GB.

### get_zip_contents(path: str, zip_relative_path: str = '', recursive: bool = False, region: str | None = None) → [ZipArchiveInfo](#inductiva.storage.storage.ZipArchiveInfo)

Retrieve the contents of a ZIP archive from a given path.

* **Parameters:**
  * **path** (*str*) – The full path to the ZIP archive.
  * **zip_relative_path** (*str* *,* *optional*) – A relative path inside the ZIP
    archive to filter the contents. Defaults to an empty string,
    which lists all files within the archive.
  * **recursive** (*bool* *,* *optional*) – If True, list contents recursively within
    the specified zip_relative_path. If False, list only top-level
    files and directories within the specified zip_relative_path.
    Defaults to False.
  * **region** (*str* *,* *optional*) – The region of the remote storage. If not
    specified, the user’s default region is assumed.
* **Returns:**
  An object containing the total size of the ZIP archive
  : and a list of ZipFileInfo objects representing the files
    within the specified ZIP archive.
* **Return type:**
  [ZipArchiveInfo](#inductiva.storage.storage.ZipArchiveInfo)

### listdir(path='/', region: str | None = None, max_results: int | None = 10, order_by: Literal['size', 'creation_time'] = 'creation_time', sort_order: Literal['asc', 'desc'] = 'desc', recursive: bool = False, print_results: bool = True)

List and display the contents of the user’s storage.
:param path: Storage directory to list. Default is root.
:type path: str
:param region: The storage region to query. If omitted, the

> user’s storage in the default region is return. Specify “all” to
> include storage from every available region (applies only to users
> with storage in multiple regions).
* **Parameters:**
  * **max_results** (*int*) – The maximum number of results to return. If not set,
    all entries are returned.
  * **order_by** (*str*) – The field to sort the contents by.
  * **sort_order** (*str*) – Whether to sort the contents in ascending or
  * **order.** (*descending*)
  * **recursive** (*bool*) – Flag to get the size and creation time of
    subdirectories.
  * **print_results** (*bool*) – Flag to print storage table.
* **Returns:**
  A list of dictionaries containing information about
  the size, the name and the creation time of each content that can
  easily be converted to a dataframe.
* **Return type:**
  list of dict

This function prints a table with the storage content information:
: Name            Size            Creation Time
  1234            5.68 MiB        29 Sep, 14:12:00
  12345           374.85 KiB      29 Sep, 14:13:10
  1234567         97.59 KiB       29 Sep, 14:13:24
  123             0 B             29 Sep, 14:13:29
  123456          0 B             29 Sep, 14:13:18

You can use this information to delete the contents you don’t need
anymore and further inspect task outputs and logs using the Task
class.

### multipart_upload(path, parts_size, upload_parts, complete_multipart_url, region: str | None = None)

Perform the multipart upload using the server.

### remove(remote_path: str, region: str | None = None)

Removes a file or directory from the remote location.

Parameters:
- remote_path (str): The path to the remote file or directory.

### upload(local_path: str, remote_dir: str, region: str | None = None)

Upload a local file or directory to a specified remote directory.

* **Parameters:**
  * **local_path** (*str*) – The path to the local file or directory to be
    uploaded.
  * **remote_dir** (*str*) – The remote directory where the file will
    be uploaded.
  * **region** (*str* *,* *optional*) – The region of the remote storage. If not
    specified, the user’s default region is assumed.

### Example

Upload a file to a remote directory:

```python
inductiva.storage.upload('local/path/file.txt', 'my_data')
```

Upload a directory to a remote location:

```python
inductiva.storage.upload('local/path/folder', 'my_data')
```

### upload_from_url(url: str, remote_dir: str, file_name: str | None = None, region: str | None = None)

Upload a file from a given URL to a specified remote directory.

If no file name is provided, the function extracts the name from the URL.

* **Parameters:**
  * **url** (*str*) – The URL of the file to be uploaded.
  * **remote_dir** (*str*) – The path to the remote directory where the file will
    be stored.
  * **file_name** (*str* *,* *optional*) – The name to save the uploaded file as.
    If not specified, the name will be extracted from the URL.
  * **region** (*str* *,* *optional*) – The region of the remote storage. If not
    specified, the user’s default region is assumed.

::docsbannersmall
::
