# Download + Delete


For a list of tasks that take up too much storage and drain credits, use a
simple loop to download the output files and then remove them from your cloud
storage.

```python
import inductiva

for task in big_tasks:
    task.download_outputs()
    task.remove_remote_files()
```

It's this simple.

> **Note**: Keep in mind that downloading the files has a cost (learn more about costs [here](https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost)). If you can, download only the needed files instead of the whole storage. (learn more how you can do that [here](../../download-file-from-project))
