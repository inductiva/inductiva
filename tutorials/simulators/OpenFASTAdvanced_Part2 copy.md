---
orphan: true
---

# Advanced OpenFAST Tutorial: Part 2

## Recap: Running a Single OpenFAST Simulation  

In the first part of this tutorial, we explored how to set up and run a single
OpenFAST simulation using Inductiva. We covered the necessary file preparations,
built the required DLL, and executed the simulation on a cost-effective cloud
machine. While running a single case on the cloud may not always be the best
option due to CPU clock speed limitations, the real advantage of Inductiva lies
in its ability to scale simulations effortlessly. Now, in Part 2, we’ll leverage
this power to run hundreds of OpenFAST simulations in parallel, drastically
reducing total computation time.

## From 1 to 100
Inductiva does not help run one OpenFast simulation faster, but
it helps you run many simulations in parallel. So, let's assume
that you need to study the impact of changing a certain 
input paramter of your simulation. For the sake of demonstration,
let us assume that we want to study what happens when the off-shore
turbine of this example is installed in locations with different 
water depths. 

In the "Enviromental Conditions" section of the `5MW_OC4Semi_WSt_WavesWN.fst`
parameter file one can see that parameter WtrDpth has been set to 200 m:


```
---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
    9.80665   Gravity         - Gravitational acceleration (m/s^2)
      1.225   AirDens         - Air density (kg/m^3)
       1025   WtrDens         - Water density (kg/m^3)
  1.464E-05   KinVisc         - Kinematic viscosity of working fluid (m^2/s)
        335   SpdSound        - Speed of sound in working fluid (m/s)
     103500   Patm            - Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
       1700   Pvap            - Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
        200   WtrDpth         - Water depth (m)
          0   MSL2SWL         - Offset between still-water level and mean sea level (m) [positive upward]
```

We are going to use Inductiva API to run variations of this
base simulation where we set the WtrDpth from 180 to 220 meter,
at 1 meter steps. This means we are going to run 40 simulations.
But we are going to run these 40 simulation in parallel.

### Step 1: Parametrize the input file `5MW_OC4Semi_WSt_WavesWN.fst`

Inductiva lets you transform fixed parameters in your simulation
configuration files into variables that you can set programmatically
via Python scripting. That is, we will be able to change the `WtrDpth`
defined in the `5MW_OC4Semi_WSt_WavesWN.fst` input file.

To do that you need to edit your `5MW_OC4Semi_WSt_WavesWN.fst` from this:

```
---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
...
        200   WtrDpth         - Water depth (m)
...
```

To this:

```
---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
...
        {{ water_depth }}   WtrDpth         - Water depth (m)
...
```

After doing this small edit you will have to save your input file with the
following name `5MW_OC4Semi_WSt_WavesWN.fst.jinja`
(watch the ".jinja" extension). This will let Inductiva's templating 
engine know that a Python variable named "water_depth" should be 
used to set the right scalar value in `5MW_OC4Semi_WSt_WavesWN.fst`.

How is this done in practice? It's very easy. The script below
shows how we can now set the value of the WtrDpth parameter from
Python, and run a variation of the original simulation for a 
water depth of 190 meters, instaed of 200 meters:

```python
import inductiva

# Allocate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="n2-highcpu-2",
    spot=True
)

water_depth = 190

print(f"Preparing files for depth = {water_depth}")
target_dir = f"variations/params_for_depth_{water_depth}"

# This is where we make the substitution and set
# the value of WtrDpth in the modified config file
# 5MW_OC4Semi_WSt_WavesWN.fst.jinja
inductiva.TemplateManager.render_dir(
    source_dir="openfast-5MW_OC4Semi_WSt_WavesWN",
    target_dir=target_dir,
    overwrite=True,
    water_depth=water_depth)

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST(
    version="4.0.2")

task = openfast.run(
    input_dir=target_dir,
    commands=[
        "openfast 5MW_OC4Semi_WSt_WavesWN/"
        "5MW_OC4Semi_WSt_WavesWN.fst"],
    on=cloud_machine)

task.wait()

cloud_machine.terminate()
```

That's it!

Now, the good thing about Inductiva API is that it is a
Python API, so you can literally just do a for loop to
iterate over all the range of values for WtrDpth. That's 
what we are going to do next.

### Step 2: Writing a "for loop" and adding more machines


#### Overview
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
    max_machines=40)

openfast_project = inductiva.projects.Project(
    name="Openfast_WavesWN",
    append=True)
openfast_project.open()

my_tasks = []
for depth in range(180, 220):

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

    my_tasks.append(task)

while(not all([x.is_terminal() for x in my_tasks])):
    print("Waiting for ALL tasks to finish")
    print(f"Finished: {sum([x.is_terminal() for x in my_tasks])}")
    time.sleep(5)

openfast_project.close()
cloud_machine.terminate()
```

Let's now break this script into parts.

#### Allocating the Cloud Machine Group  

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
    num_machines=40)
```

#### Setting Up Our Project  

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

#### Preparing and Running the Simulation  

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
for depth in range(180, 220):

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

    # Store the task for later synchronization
    my_tasks.append(task)
```

#### Waiting for the simulations to finish

In this last part of the code we are just waiting for all the tasks to finish.

Once all tasks are done we close our project and turn off our cloud machines.

```python
while(not all([x.is_terminal() for x in my_tasks])):
    print("Waiting for ALL tasks to finish")
    print(f"Finished: {sum([x.is_terminal() for x in my_tasks])}")
    time.sleep(5)

openfast_project.close()
cloud_machine.terminate()
```

>Note: You can list all the tasks from this project by running this command on
your command line interface `inductiva tasks list -p Openfast_WavesWN -n 40`

```
Showing tasks for project: Openfast_WavesWN.

 ID                          SIMULATOR     STATUS     SUBMITTED         STARTED           COMPUTATION TIME     RESOURCE TYPE
 9ousr15npia4f9fphfjm1rmrn   openfast      Success    12/02, 15:25:00   12/02, 15:25:00   0:00:34              GCP n2-highcpu-2
 pt3i70i4hmzn740og3ber5yvs   openfast      Success    12/02, 15:24:54   12/02, 15:24:54   0:00:32              GCP n2-highcpu-2
 60ygm64qsd56fmcjxgytq93pr   openfast      Success    12/02, 15:24:47   12/02, 15:24:47   0:00:35              GCP n2-highcpu-2
 sy2a5r3cby0in6151ze06y29g   openfast      Success    12/02, 15:24:39   12/02, 15:24:39   0:00:35              GCP n2-highcpu-2
 snj7xr01k099qbh9ytf97qtbd   openfast      Success    12/02, 15:24:32   12/02, 15:24:32   0:00:37              GCP n2-highcpu-2
 h4qtfsq43si18a6oyvj74qrce   openfast      Success    12/02, 15:24:27   12/02, 15:24:27   0:00:35              GCP n2-highcpu-2
 j85bg918c1i0obj1izumtwwzw   openfast      Success    12/02, 15:24:19   12/02, 15:24:19   0:00:32              GCP n2-highcpu-2
 jxn9a281jkyi4f62n4w80a1yp   openfast      Success    12/02, 15:24:15   12/02, 15:24:15   0:00:31              GCP n2-highcpu-2
 tyqb7rnsyzi94wqgcn58cc2y6   openfast      Success    12/02, 15:24:08   12/02, 15:24:08   0:00:35              GCP n2-highcpu-2
 1d61519z29drg2jp4nhknauoo   openfast      Success    12/02, 15:24:00   12/02, 15:24:00   0:00:30              GCP n2-highcpu-2
 68ksaoricizknyhyg65r6cuqo   openfast      Success    12/02, 15:23:55   12/02, 15:23:55   0:00:35              GCP n2-highcpu-2
 3qoor8iooxtmmm5dozs0o40px   openfast      Success    12/02, 15:23:48   12/02, 15:23:48   0:00:32              GCP n2-highcpu-2
 fr286c6cfrmcus12nh4n37m78   openfast      Success    12/02, 15:23:42   12/02, 15:23:42   0:00:34              GCP n2-highcpu-2
 b0mhrbzy9fwrfuj3f84q8ou75   openfast      Success    12/02, 15:23:35   12/02, 15:23:35   0:00:32              GCP n2-highcpu-2
 n0jb7r2rpuyik5xr7e96d024r   openfast      Success    12/02, 15:23:27   12/02, 15:23:27   0:00:36              GCP n2-highcpu-2
 d6ght4mi7u10mzkstm33fpvty   openfast      Success    12/02, 15:23:20   12/02, 15:23:20   0:00:33              GCP n2-highcpu-2
 pqv4qmy2so820yoqm6dj51na0   openfast      Success    12/02, 15:23:13   12/02, 15:23:13   0:00:30              GCP n2-highcpu-2
 eyw55mcvecqhjxnwncpy6wqkm   openfast      Success    12/02, 15:23:05   12/02, 15:23:05   0:00:35              GCP n2-highcpu-2
 t98ihj29ieep26ndstns63u9r   openfast      Success    12/02, 15:22:57   12/02, 15:22:57   0:00:35              GCP n2-highcpu-2
 qgehpoemdijghnu27hqjed813   openfast      Success    12/02, 15:22:49   12/02, 15:22:49   0:00:31              GCP n2-highcpu-2
 bfn4ok5cwy93bmi1hg6vvk0rs   openfast      Success    12/02, 15:22:44   12/02, 15:22:44   0:00:32              GCP n2-highcpu-2
 4jho5mrt1dvr1c4wcsq9hje2z   openfast      Success    12/02, 15:22:37   12/02, 15:22:37   0:00:35              GCP n2-highcpu-2
 r1kqwy76fvus9zaevv9xzzc02   openfast      Success    12/02, 15:22:30   12/02, 15:22:30   0:00:33              GCP n2-highcpu-2
 tif6v2qugbkcjcum2hna5ttow   openfast      Success    12/02, 15:22:21   12/02, 15:22:21   0:00:37              GCP n2-highcpu-2
 smnphh002n0pxqnte7r6ihekg   openfast      Success    12/02, 15:22:14   12/02, 15:22:14   0:00:33              GCP n2-highcpu-2
 hmf133908834ms1hg9oh49q4e   openfast      Success    12/02, 15:22:06   12/02, 15:22:06   0:00:33              GCP n2-highcpu-2
 k6xoup4jqczhjuwmiig4c7eo9   openfast      Success    12/02, 15:22:02   12/02, 15:22:02   0:00:33              GCP n2-highcpu-2
 u3my9tc6deh72xan9j7kn0syj   openfast      Success    12/02, 15:21:55   12/02, 15:21:55   0:00:35              GCP n2-highcpu-2
 4rzzrir5dmhu4g49ft2t8dmw0   openfast      Success    12/02, 15:21:46   12/02, 15:21:46   0:00:32              GCP n2-highcpu-2
 7jkebn05rfsfcprto45k8nb4h   openfast      Success    12/02, 15:21:39   12/02, 15:21:39   0:00:33              GCP n2-highcpu-2
 4vocm7rjiehyiyb05mmnpxw0u   openfast      Success    12/02, 15:21:31   12/02, 15:21:31   0:00:33              GCP n2-highcpu-2
 kndw77s3jzumypsnssszp4m4h   openfast      Success    12/02, 15:21:23   12/02, 15:21:23   0:00:34              GCP n2-highcpu-2
 btctwrbp2z69db7w5rbdyuatm   openfast      Success    12/02, 15:21:17   12/02, 15:21:17   0:00:33              GCP n2-highcpu-2
 bzhd49wxjkroc30f4vgs8u4o8   openfast      Success    12/02, 15:21:09   12/02, 15:21:09   0:00:33              GCP n2-highcpu-2
 qfrv5vkims0zqbv7gevviiofg   openfast      Success    12/02, 15:21:02   12/02, 15:21:02   0:00:30              GCP n2-highcpu-2
 691acms5avqwvpcqf8kvecw33   openfast      Success    12/02, 15:20:54   12/02, 15:20:54   0:00:35              GCP n2-highcpu-2
 5m9nnjw28h31puxybbzk1mq6h   openfast      Success    12/02, 15:20:47   12/02, 15:20:47   0:00:31              GCP n2-highcpu-2
 8gvt5i0x0qeaj2x7ftf7pwxyd   openfast      Success    12/02, 15:20:37   12/02, 15:20:37   0:00:31              GCP n2-highcpu-2
 gnq3smc96ht404hz4ucvoovw6   openfast      Success    12/02, 15:20:14   12/02, 15:20:19   0:00:33              GCP n2-highcpu-2
 iox5b1040tagfnwi9a3j26axr   openfast      Success    12/02, 15:20:07   12/02, 15:20:19   0:00:31              GCP n2-highcpu-2

To see more details about a task, use `inductiva tasks info <task_id>`.
```

## Conclusion  

In this tutorial, we demonstrated how to leverage Inductiva to efficiently run
multiple OpenFAST simulations in parallel. While a single OpenFAST simulation
may run faster on a high-frequency local machine, scaling up to dozens or
hundreds of simulations is where Inductiva truly shines. By automating input
file modifications and distributing simulations across inexpensive cloud
machines, you can drastically reduce total computation time and streamline
large-scale studies.  

With just a few lines of Python, you can scale your OpenFAST projects
effortlessly—saving time, optimizing resources, and accelerating your research. 🚀