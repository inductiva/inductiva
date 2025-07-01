# Setting a TTL (Time-to-Live) on Your Simulations

To set a TTL (Time-to-Live) for your simulation, you need to use the `time_to_live` 
parameter when submitting the task. This parameter defines the maximum duration the 
task is allowed to run, expressed as a string with a time unit suffix. Supported formats 
include minutes (e.g., `"10m"`) or hours (e.g., `"2h"`). Once the specified time has 
elapsed from the task's start, the simulation will be automatically terminated. This 
feature helps you control resource usage and avoid running tasks longer than necessary 
(e.g., perfect for quick test runs).

```python
import inductiva

# ...

task = simulator.run(
    input_dir=input_dir,
    sim_config_filename=sim_config_filename,
    on=cloud_machine,
    time_to_live="1m",  # Sets TTL to 1 minute
)
```


You'll see logs like these in your terminal:

```bash
■ Your task exceeded its configured time-to-live (TTL) and was automatically stopped.
Downloading stdout and stderr files to inductiva_output/h04n8swuji8jgrex5sq7smsg2/outputs...
Partial download completed to inductiva_output/h04n8swuji8jgrex5sq7smsg2/outputs.
```

In the [Inductiva Console](https://console.inductiva.ai/dashboard), open the 
[Tasks](https://console.inductiva.ai/tasks) tab and click on your task to view its details. 
You'll see something like this:

![TTL task details in Console](static/console-ttl.png)
