# Running multiple simulations in parallel

**Inductiva API** allows you to launch multiple simulations in parallel. After a simulation is submitted, users are free to execute other commands, since the simulation runs independently in the background. Users are not limited by waiting for the simulation to finish.

Users can prepare several simulations and launch them all together in parallel. Let's look at an example using the wind tunnel scenario:

```python
import inductiva

vehicle_url = "https://raw.githubusercontent.com/inductiva/inductiva/main" \
              "/assets/vehicle.obj"
vehicle_path = inductiva.utils.files.download_from_url(vehicle_url)

task_list = []
velocity_list=[1, 10, 20, 30, 40]
for velocity in velocity_list:
  scenario = inductiva.fluids.WindTunnel(flow_velocity=[velocity, 0, 0])
  task = scenario.simulate(object_path=vehicle_path)
  task_list.append(task)
```

All simulations are launched in one go, allowing users to continue working on other things. To monitor the progress of individual simulations, users can use `task.get_status()`, or they can view a list of the most recent tasks launched by using `inductiva.tasks.list(last_n=5)`.

Finally, to retrieve the results the user can use `task.download_outputs()`, which waits for the simulation to finish before downloading the results. During this waiting period, it temporarily blocks the execution of other code. Check the [Tasks section](https://github.com/inductiva/inductiva/tree/main/inductiva/tasks#tasks) for more information on how to do this.