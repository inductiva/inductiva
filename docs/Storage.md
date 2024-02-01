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
This will enumerate the contents of the root directory and produce a table showcasing the storage directories: 

```bash

       NAME                       SIZE            CREATION TIME
       1699461562982775346/       711.82 MB       08 Nov, 16:39:23
       1699461561456065056/       706.74 MB       08 Nov, 16:39:21
       1698772828190838710/       706.44 MB       31 Oct, 17:20:28
       1698772827312137069/       706.18 MB       31 Oct, 17:20:27
       1699461562181988760/       699.30 MB       08 Nov, 16:39:22
       1698772829464378916/       698.80 MB       31 Oct, 17:20:29
       1699461560739645733/       693.81 MB       08 Nov, 16:39:20
       1698772828793721336/       688.82 MB       31 Oct, 17:20:28
       1698751109500351563/       685.33 MB       31 Oct, 11:18:29
       1699461560029476871/       673.81 MB       08 Nov, 16:39:20
```
To examine the contents of a specific folder, execute: 

```python
import inductiva
inductiva.storage.listdir(path = "1234", max_results=10, order_by="size", sort_order="desc")
```

The `order_by` argument allows you to sort the table by size or creation date, while the `sort_order` argument determines whether the list is displayed in ascending or descending order. 

### Removing directories

The table provides valuable information that can guide your decision to remove certain directories. This can be accomplished using the `inductiva.storage.rmdir` function. 

```python
import inductiva
inductiva.storage.rmdir(path="1234")
```
Running the above code will permanently delete the directory named "1234" from your remote storage. Please note, this action is irreversible and the deleted directory cannot be restored.

