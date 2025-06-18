# Set up an Elastic Machine Group

The `ElasticMachineGroup`, similarly to the standard `MachineGroup`, is composed of 
a group of homogeneous machines that work individually to run multiple simulations. 
The difference is that the number of active machines is scaled up and/or down 
automatically based on the simulations in the queue. This prevents computational 
resources from being idle when there are no sufficient simulations to run and 
allows scaling up the computational resources when the queue is full.

Note that, the elasticity is independent of each machine being preemptible, i.e.,
these can be combined and Inductiva API manages the simulations accordingly.

Users can configure the resource with the minimum number of active machines they
want to always be active, `min_machines`, and the maximum number of machines that
can be active at all times, `max_machines`.

Let's we start an `ElasticMachineGroup` object:

```python
import inductiva

# Configure an elastic machine group to start with a minimum of 1 machine up to a
# maximum of 3, each with a data disk of 30 Gb.
elastic_machine_group = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="c2-standard-30",
    min_machines=1,
    max_machines=3,
    data_disk_gb=30)

# Launch the Elastic machine group to make it available to run simulations:
elastic_machine_group.start()
```

Once started, simulations can be submitted to the queue of the elastic machine group.

To explore our elastic machine group, we will follow the example of 
[running multiple simulations in parallel](./run-parallel_simulations.md), based on the 
[templating mechanism](https://tutorials.inductiva.ai/intro_to_api/templating.html)
built in the Inductiva API, but now with a scalable infrastructure.

```python
import inductiva

# Download the input files for the SWASH simulation
template_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-template-example.zip", unzip=True)

# Initialize the SWASH simulator
swash = inductiva.simulators.SWASH()

# Explore the simulation for different water levels
water_levels_list = [3.5, 3.75, 4.0, 4.5, 5.0]

# Launch multiple simulations
for i, water_level in enumerate(water_levels_list):
    target_dir = f"./inductiva_input/swash-sim-{i}"  
    inductiva.TemplateManager.render_dir(
                            source_dir=template_dir,
                            target_dir=target_dir,
                            water_level=water_level,
                            overwrite=False)

    # Run the simulation on the dedicated MachineGroup
    task = swash.run(
        input_dir=target_dir,
        sim_config_filename="input.sws",
        on=elastic_machine_group)
```

As our simulations are submitted to the queue of the elastic machine group, we
will use the CLI to check the status of our simulations and the scaling of the resources.

> **__TIP__ 1** Open two terminal windows to monitor both the simulations and resources
simultaneously.

> **__TIP__ 2** With the watch method available in some OS's, one can repeatedly apply
the CLI commands to monitor the tasks and the resources.

To track the active resources we use
```bash
$ watch inductiva resources ls
```
and for our last 5 tasks submitted, we tracked them with
```bash
$ watch inductiva tasks list -n 5
```

Right after launch, the elastic machine group starts with a single active machine
and picks up one simulation:

**Resources monitoring**
```bash
Every 2.0s: inductiva resources ls

Active Resources:

       NAME                                MACHINE TYPE         ELASTIC         TYPE           # MACHINES         DATA SIZE IN GB         SPOT         STARTED AT (UTC)
       api-t83rivz4yeh29k3hy4dofa38d       c2-standard-30       True            standard       1/3                30                      False        07 Feb, 11:47:20
```

**Tasks monitoring**
```bash
Every 2.0s: inductiva tasks list                                                                                                                   
            
       ID                              SIMULATOR         STATUS          SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       r4kerxf4b53krgn0s3fyece3b       swash             submitted       07 Feb, 11:47:49       n/a                    n/a                      n/a
       j9qzrpiohgt7x97od3tw4wccd       swash             submitted       07 Feb, 11:47:48       n/a                    n/a                      n/a
       iqi71gonoacfj7fknox3rvnq2       swash             submitted       07 Feb, 11:47:46       n/a                    n/a                      n/a
       dxmnxdrfrv84pfbzbvm9v0dat       swash             submitted       07 Feb, 11:47:44       n/a                    n/a                      n/a
       bgtwgnnyq5qa5hecegzdx6okr       swash             started         07 Feb, 11:47:42       07 Feb, 11:48:08       *0:01:50                 c2-standard-30
```

Notice that on the tasks, there are 4 simulations submitted and
waiting on the queue of the elastic machine group. With the awareness of a non-empty
queue, after 2 minutes another machine becomes active and after 4 minutes, the
elastic machine group is fully active and three tasks are running simultaneously:

**Resources monitoring**
```bash
Active Resources:

       NAME                                MACHINE TYPE         ELASTIC         TYPE           # MACHINES         DATA SIZE IN GB         SPOT         STARTED AT (UTC)
       api-t83rivz4yeh29k3hy4dofa38d       c2-standard-30       True            standard       3/3                30                      False        07 Feb, 11:47:20
```

**Tasks monitoring**
```bash
Every 2.0s: inductiva tasks list                                                                                                     Ivans-MacBook-Air.local: Wed Feb  7 11:54:12 2024


       ID                              SIMULATOR         STATUS          SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       r4kerxf4b53krgn0s3fyece3b       swash             submitted       07 Feb, 11:47:49       n/a                    n/a                      n/a
       j9qzrpiohgt7x97od3tw4wccd       swash             submitted       07 Feb, 11:47:48       n/a                    n/a                      n/a
       iqi71gonoacfj7fknox3rvnq2       swash             started         07 Feb, 11:47:46       07 Feb, 11:52:44       *0:01:31                 c2-standard-30
       dxmnxdrfrv84pfbzbvm9v0dat       swash             started         07 Feb, 11:47:44       07 Feb, 11:50:27       *0:03:50                 c2-standard-30
       bgtwgnnyq5qa5hecegzdx6okr       swash             started         07 Feb, 11:47:42       07 Feb, 11:48:08       *0:06:10                 c2-standard-30
```

Now, as the tasks finish running, the still active machine pick-up another task from
the queue until all are completed. When all are complete the elastic machine group will
start scaling down until it stays with the `min_machines`, in our case, 1 machine.

**Tasks monitoring**
```bash
Every 2.0s: inductiva tasks list                                                                                        

       ID                              SIMULATOR         STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       r4kerxf4b53krgn0s3fyece3b       swash             success        07 Feb, 11:47:49       07 Feb, 12:00:55       0:10:29                  c2-standard-30
       j9qzrpiohgt7x97od3tw4wccd       swash             success        07 Feb, 11:47:48       07 Feb, 11:58:10       0:10:03                  c2-standard-30
       iqi71gonoacfj7fknox3rvnq2       swash             success        07 Feb, 11:47:46       07 Feb, 11:52:44       0:10:02                  c2-standard-30
       dxmnxdrfrv84pfbzbvm9v0dat       swash             success        07 Feb, 11:47:44       07 Feb, 11:50:27       0:10:20                  c2-standard-30
       bgtwgnnyq5qa5hecegzdx6okr       swash             success        07 Feb, 11:47:42       07 Feb, 11:48:08       0:09:54                  c2-standard-30
```

**Resources monitoring**
```bash
Every 2.0s: inductiva resources ls                                                                                         

Active Resources:

       NAME                                MACHINE TYPE         ELASTIC         TYPE           # MACHINES         DATA SIZE IN GB         SPOT         STARTED AT (UTC)
       api-t83rivz4yeh29k3hy4dofa38d       c2-standard-30       True            standard       3/3                30                      False        07 Feb, 11:47:20
```

The elastic machine group plays the trade-off between the idle time of resources and 
the queue time of the simulations without extra configuration by the user. 

At the end, there will always be the number of `min_machines` active,
therefore, don't forget to terminate your resources with:
```bash
$ inductiva resources terminate --all -y
```
