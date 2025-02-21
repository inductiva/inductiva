---
orphan: true
---

# Running 50 Simulations

In the previous tutorial, we set up a templating system that allows us to
modify the `WtrDpth` parameter in the OpenFAST input file programmatically.
This enables us to, easly, generate multiple simulation configurations with different
water depths.

Now, we will take the next step and use a for loop to automate the process,
launching 50 simulations in parallel. This approach demonstrates the true power
of Inductiva: scaling simulations efficiently to save computation time.

Let's dive in!

## Writing a "for loop" and adding more machines


### Overview
Here is the code to run all the simulations in parallel.

```python
import inductiva
import time

# Allocate cloud machine
cloud_machine = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="n2-highcpu-2",
    spot=True,
    min_machines=1,
    max_machines=50)

openfast_project = inductiva.projects.Project(
    name="Openfast_WavesWN",
    append=True)
openfast_project.open()

for depth in range(100, 200, 2):

    print(f"Preparing files for depth = {depth}")
    target_dir = f"variations/params_for_depth_{depth}"

    inductiva.TemplateManager.render_dir(
        source_dir="input_files",
        target_dir=target_dir,
        overwrite=True,
        water_depth=depth)

    # Initialize OpenFAST simulator
    openfast = inductiva.simulators.OpenFAST(
        version="4.0.2",
        use_dev=True)

    task = openfast.run(
        input_dir=target_dir,
        commands=[
        "openfast 5MW_OC4Semi_WSt_WavesWN/"
        "5MW_OC4Semi_WSt_WavesWN.fst"],
        on=cloud_machine)

openfast_project.wait()
openfast_project.close()
cloud_machine.terminate()
```

Let's now break this script into parts.

### Allocating the Cloud Machine Group  

We are allocating an elastic machine group to run our simulations. This group
has a minimum number of machines active and a maximum capacity. Machines are
dynamically turned on or off based on demand, ensuring efficient resource utilization.  

Since some simulations may take longer than others, this elastic setup prevents
idle waiting. When a machine completes its task and has no more simulations to
run, it does not need to wait for all other machines to finish. Instead, the
elastic machine group automatically scales down, optimizing cloud usage and
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

### Setting Up Our Project  

To keep our simulations organized, we use a project-based approach. Setting up
a project is straightforward: we create it with **append mode enabled**, which
indicates that we’ll be adding tasks to the project rather than just querying
its data. Once the project is created and opened, all subsequent simulations
will automatically be associated with our `Openfast_WavesWN` project.  

```python
openfast_project = inductiva.projects.Project(
    name="Openfast_WavesWN",
    append=True
)
openfast_project.open()
```  

This ensures that all simulations remain structured and easily accessible
within the project.

### Preparing and Running the Simulation  

We now enter a loop responsible for generating new input files and running each
simulation.  

First, we iterate over the specified water depths. For each depth, we generate
input files using the `render_dir` method, which replaces the `water_depth`
variable in our `5MW_OC4Semi_WSt_WavesWN.fst.jinja` template with the
corresponding depth value. The generated input files are then stored in a
dedicated folder, `target_dir` (`variations/params_for_depth_{depth}`).  

Next, we initialize the OpenFAST simulator and run the simulation using the
newly created input directory. Each simulation task is stored in a list,
allowing us to track and wait for all tasks to complete.  

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
        on=cloud_machine
    )
```

### Waiting for the simulations to finish

In this last part of the code we are just waiting for all the tasks to finish.

Once all tasks are done we close our project and turn off our cloud machines.

```python
openfast_project.wait()

openfast_project.close()
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

## Conclusion  

Previously, we demonstrated how to leverage Inductiva to efficiently run
multiple OpenFAST simulations in parallel. While a single OpenFAST simulation
may run faster on a high-frequency local machine, scaling up to dozens or
hundreds of simulations is where Inductiva truly shines. By automating input
file modifications and distributing simulations across inexpensive cloud
machines, you can drastically reduce total computation time and streamline
large-scale studies.  

With just a few lines of Python, you can scale your OpenFAST projects
effortlessly—saving time, optimizing resources, and accelerating your research. 🚀

[Downloading The Results](OpenFASTAdvanced_Part6.md)
