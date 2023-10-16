# Storage management

The Inductiva API provides you with tools to manage your remote storage effectively. 
With the Inductiva storage module, you can easily navigate your storage, evaluate the space used, and delete specific directories as needed.
You have the ability to organize your storage by specifying the directory where your simulation outputs should be saved. The directory containing these outputs is automatically named after the task ID. Check [tasks](https://github.com/inductiva/inductiva/tree/main/inductiva/tasks) for more information about this.
Let's illustrate this with some examples:

## Usage example
### Determining the amount of storage in use

```python
import inductiva
space_used = inductiva.storage.get_space_used()
```
The above code snippet will return the total amount of storage space currently being utilized in your storage.

### Viewing storage contents
After determining the total storage space used, you may want to identify which directories are consuming the most storage. This can be achieved using the `inductiva.storage.listdir` function:

```python
import inductiva
inductiva.storage.listdir(max_results=10, order_by="size", sort_order="desc")
```
This function generates a table displaying the storage contents: 

```markdown
| Name   | Size       | Creation Time   |
|--------|------------|-----------------|
| 1234   | 5.68 MiB   | 29 Sep, 14:12:00|
| 12345  | 374.85 KiB | 29 Sep, 14:13:10|
| 1234567| 97.59 KiB  | 29 Sep, 14:13:24|
| 123    | 0 B        | 29 Sep, 14:13:29|
```

The `order_by` argument allows you to sort the table by size or creation date, while the `sort_order` argument determines whether the list is displayed in ascending or descending order. 

### Removing directories

The table provides valuable information that can guide your decision to remove certain directories. This can be accomplished using the `inductiva.storage.rmdir` function. 

```python
import inductiva
inductiva.storage.rmdir(path="1234")
```
Running the above code will permanently delete the directory named "1234" from your remote storage. Please note, this action is irreversible and the deleted directory cannot be restored.
