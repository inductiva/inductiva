# Storage Management

The Inductiva API provides you with the tools to manage your cloud storage effectively. Currently, we support Google Cloud, but we're working on expanding our services to include multiple cloud platforms soon!
With the Inductiva storage module, you can easily navigate your storage bucket, evaluate the space used, and delete specific directories as needed.
Your storage bucket is organized by task ID, with each task corresponding to a folder named after its ID. These folders contain both the input and output files for each simulation.
These functions are closely linked with the task class, allowing for seamless manipulation of task outputs.
Let's illustrate this with some examples:

## Usage Example
### Determining the Amount of Storage Used

```python
import inductiva
space_used = inductiva.storage.get_space_used()
```
The above code snippet will return the total amount of storage space currently being utilized in your bucket.

### Viewing Storage Contents
After determining the total storage space used, you may want to identify which directories are consuming the most storage. This can be achieved using the `inductiva.storage.listdir` function:

```python
import inductiva
inductiva.storage.listdir(max_results=10, order_by="size", sort_order="desc")
```
This function generates a table displaying the storage content information:
        Name            Size            Creation Time
        1234            5.68 MiB        29 Sep, 14:12:00
        12345           374.85 KiB      29 Sep, 14:13:10
        1234567         97.59 KiB       29 Sep, 14:13:24
        123             0 B             29 Sep, 14:13:29
        
The `order_by` argument allows you to sort the table by size or creation date, while the `sort_order` argument determines whether the list is displayed in ascending or descending order. 

### Removing Directories

The table provides valuable information that can guide your decision to remove certain directories. This can be accomplished using the `inductiva.storage.rmdir` function. 

```python
import inductiva
inductiva.storage.rmdir(path="1234")
```
Executing the above code will effectively remove the directory named "1234" from your remote storage. 
