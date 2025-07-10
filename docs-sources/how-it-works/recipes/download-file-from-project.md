# How to Download a Specific File from All Tasks in a Project  

## The Challenge  

You've just ran hundreds of simulations, each generating multiple output files.
Now, you need to collect a specific file from every simulation. Manually
searching, downloading, and organizing these files is not only tedious but also
prone to errors.  

Instead of doing this manually, why not automate the process?

With the right approach, you can gather these files efficiently, ensuring
accuracy while saving valuable time.  

## The Solution  

Inductiva makes this process simple. The following script demonstrates how to
automatically download a specific file from every task in your project:  

```python
import inductiva

project = inductiva.projects.Project(
    name="my-project")

for task in project.get_tasks():
    # Download a specific output file from each task
    task.download_outputs(["path/to/output/file.txt"])
```

You can also download multiple files at once by providing a list of file paths.
This way, managing large-scale simulations becomes effortless and error-free.