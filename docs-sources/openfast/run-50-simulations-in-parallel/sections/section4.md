# Run 50 Simulations
In the previous part of this tutorial we set up a templating system that allows us to programmatically modify the `WtrDpth` parameter in the OpenFAST input file. 
This allows us to easily generate multiple simulation configurations with different water depths.

We will now take the next step and use a for loop to automate the process and run 50 simulations in parallel. This approach demonstrates the true power 
of Inductiva: efficiently scaling simulations to save computation time.

## Writing a “for loop” and adding more machines
We will first show you all the necessary code and then analyze it step by step. Each section of code is easy to understand and builds on the examples shown 
in the previous parts of this tutorial.

```python
import inductiva

# Allocate 50 cloud machines
# Note that we are using an ElasticMachineGroup
cloud_machine = inductiva.resources.ElasticMachineGroup(
   provider="GCP",
   machine_type="n2-highcpu-2",
   spot=True,
   min_machines=1,
   max_machines=50)


for depth in range(100, 200, 2):


   print(f"Preparing files for depth = {depth}")
   target_dir = f"variations/params_for_depth_{depth}"


   # This is where we instantiate the configuration
   # files with the specific values
   inductiva.TemplateManager.render_dir(
       source_dir="input_files",
       target_dir=target_dir,
       overwrite=True,
       water_depth=depth)


   # Initialize OpenFAST simulator
   openfast = inductiva.simulators.OpenFAST(
       version="4.0.2")


   task = openfast.run(
       input_dir=target_dir,
       commands=[
       "openfast 5MW_OC4Semi_WSt_WavesWN/"
       "5MW_OC4Semi_WSt_WavesWN.fst"],
       on=cloud_machine,
       project="Openfast_WavesWN")


inductiva.projects.Project("Openfast_WavesWN").wait()
cloud_machine.terminate()
```

Let's break this script into parts.

## Code Section 1: Allocating the Cloud Machine Group
We will allocate an elastic machine group to run our simulations. This group has a minimum number of active machines and a maximum capacity. 
Machines are dynamically turned on or off as demanded, ensuring efficient resource utilization.

As some simulations may take longer than others, this elastic setup prevents idle waiting. When a machine completes its task and has no more simulations to run, 
it does not have to wait for all the other machines to finish. Instead, the elastic group of machines automatically scales down, optimizing cloud usage and 
reducing costs.

```python
# Allocate cloud machine
cloud_machine = inductiva.resources.ElasticMachineGroup(
   provider="GCP",
   machine_type="n2-highcpu-2",
   spot=True,
   min_machines=1,
   max_machines=50)
```

## Code Section 2: Loop over parameter assignment and Start Simulation
We now enter a loop that is responsible for generating new input files and running each simulation.

First, we iterate over the specified water depths. For each depth, we generate input files using the `render_dir` method, which replaces the `water_depth` 
variable in our `5MW_OC4Semi_WSt_WavesWN.fst.jinja` template with the corresponding depth value. The generated input files are then stored 
in a special folder, `target_dir` (`variations/params_for_depth_{depth}`). For more details on managing template files, check out the `TemplateManager` [documentation](https://docs.inductiva.ai/en/latest/intro_to_api/templating.html).

Next, we initialize the OpenFAST simulator and run the simulation using the newly created input directory. Each simulation task is added to the project, 
allowing us to track and wait for all tasks to be completed.

```python
for depth in range(100, 200, 2):


   print(f"Preparing files for depth = {depth}")
   target_dir = f"variations/params_for_depth_{depth}"


   # Generate input files with the specified water depth
   inductiva.TemplateManager.render_dir(
       source_dir="input_files",
       target_dir=target_dir,
       overwrite=True,
       water_depth=depth
   )


   # Initialize OpenFAST simulator
   openfast = inductiva.simulators.OpenFAST(
       version="4.0.2"
   )


   # Run the simulation on the allocated cloud machines
   task = openfast.run(
       input_dir=target_dir,
       commands=[
           "openfast 5MW_OC4Semi_WSt_WavesWN/"
           "5MW_OC4Semi_WSt_WavesWN.fst"
       ],
       on=cloud_machine,
       project="Openfast_WavesWN"
   )
```

## Code Section 4: Waiting for the Simulations to finish
In this last part of the code we just wait for all the tasks to finish.

Once all the tasks are done, we turn off our cloud machines.

```python
inductiva.projects.Project("Openfast_WavesWN").wait()

cloud_machine.terminate()
```

>Note: You can list all the tasks from this project by running this command on
your command line interface `inductiva tasks list -p Openfast_WavesWN -n 50`

```
Showing tasks for project: Openfast_WavesWN.

 ID                          SIMULATOR     STATUS     SUBMITTED         STARTED           COMPUTATION TIME     RESOURCE TYPE
 9ousr15npia4f9fphfjm1rmrn   openfast      Success    12/02, 15:25:00   12/02, 15:25:00   0:00:34              GCP n2-highcpu-2
 pt3i70i4hmzn740og3ber5yvs   openfast      Success    12/02, 15:24:54   12/02, 15:24:54   0:00:32              GCP n2-highcpu-2
 60ygm64qsd56fmcjxgytq93pr   openfast      Success    12/02, 15:24:47   12/02, 15:24:47   0:00:35              GCP n2-highcpu-2
 ...
 8gvt5i0x0qeaj2x7ftf7pwxyd   openfast      Success    12/02, 15:20:37   12/02, 15:20:37   0:00:31              GCP n2-highcpu-2
 gnq3smc96ht404hz4ucvoovw6   openfast      Success    12/02, 15:20:14   12/02, 15:20:19   0:00:33              GCP n2-highcpu-2
 iox5b1040tagfnwi9a3j26axr   openfast      Success    12/02, 15:20:07   12/02, 15:20:19   0:00:31              GCP n2-highcpu-2

To see more details about a task, use `inductiva tasks info <task_id>`.
```

## Summary
We demonstrated how Inductiva can be used to efficiently run multiple OpenFAST simulations in parallel. While a single OpenFAST simulation may run faster 
on a local machine with high processing power, scaling up to dozens or hundreds of simulations is where Inductiva excels. By automating input file modifications
and distributing simulations across inexpensive cloud machines, you can drastically reduce overall computation time and streamline large-scale studies.

With just a few lines of Python, you can effortlessly scale your OpenFAST projects - saving time, optimising resources and accelerating your research. 🚀
