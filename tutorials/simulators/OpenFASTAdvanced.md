In this tutorial, we are going to show you how you can
use the Inductiva API to accelerate your OpenFast projects.

## The curious case of OpenFast

OpenFast is a pretty peculiar piece of software because
it runs only on a single thread. That means that OpenFast
does not take full advantage of the extra compute capacity
offered by modern CPUs, which support dozens, or sometimes 
hundreds, of computational threads in parallel. In fact,
the most decisive factor impacting the performance of 
OpenFast simulations is CPU clock frequency. 

Now, what is interesting about that is that, tipically, your 
500\$ home desktop has a much higher clock frequency than the 
20000\$ monster machines available in cloud centers. Your local
desktop has CPU with clock frequencies between 4 and 5GHz (but
can be overclocked to more than 5.5GHz), while cloud machines
have CPUs that will be running at around 3GHz (or less!).

But why is that? Have you noticed how fast the ventilation 
fans of your desktop PC spin when you are running some long 
simulation? Now imagine how many fans and sophisticated
cooling systems you would need if you packed 100 CPUs like
yours in a volume the size of your kitchen fridge. Hardware
density in datacenter cabinets is so high that if all the
CPUs were running at speeds of 5GHz or more, the datacenter
would be impossible to cool down.

So, if that is the case, why would you even bother using
Inductiva to run your OpenFast simulation? What advantage
could you possibly have if you decied to send your OpenFast
simulation to a cloud machine using Inductiva? 

And the answer is: **hundreds of advantages!**

Let me explain! If you are running only 1 OpenFast simulation, 
then you should run it on your desktop machine: it will be faster
there due to much higher clock speeds.  But if you need to run
hundreds or thousands of OpenFast simulations, then you can use
Inductiva to spin up hundreds of very cheap cloud machines to 
run all those simulations in parallel instead of having to run
them sequentially in your machine! This will be hundreds of times
faster. And it is super easy (and cost-effective) with Inductiva

In this tutorial, we will show you how to do this using the 
"5MW_OC4Semi_WSt_WavesWN" example,vavailable from the OpenFast 
[GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.
 
This example is an extension of the reference case described in 
["Definition of a 5-MW Reference Wind Turbine for Offshore
System Development"](https://www.nrel.gov/docs/fy09osti/38060.pdf).

All files needed are available from OpenFast GitHub so
let's get started.

## Running the baseline use case

### Step 1: Downloading you simulation files

We are now going to run the `5MW_OC4Semi_WSt_WavesWN` case as
it is originally defined in the [GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.

We start by downloading the folders `5MW_OC4Semi_WSt_WavesWN` and the `5MW_Baseline` to our local input directory. This folders are located [here](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast).

After downloading the folders and moving them your `input_files` should look something like this:

- input_files
    - 5MW_Baseline
    - 5MW_OC4Semi_WSt_WavesWN

### Step 2: Building `DISCON_OC3Hywind.dll`

After downloading the files we need to build the `DISCON_OC3Hywind.dll` file and place it in `5MW_Baseline/ServoData/`. 

To do so we need to do the following steps:
```
cd input_files/5MW_Baseline/ServoData/DISCON_OC3/
mkdir build
cd build
cmake ..
make
mv DISCON_OC3Hywind.dll ../..
```

We should now have our input directory looking like this:
  
- input_files
    - 5MW_Baseline
        - ServoData
            - DISCON_OC3Hywind.dll
    - 5MW_OC4Semi_WSt_WavesWN

> Note: you can also download de dll file [here](https://storage.googleapis.com/inductiva-simulators-sources/DISCON_OC3Hywind.dll). 
Dont forget to paste this file in the 5MW_Baseline/ServoData folder.

You now have all the necessary files to run your simulation.

## Running your simulation

For you to run this simulation you can execute the following python code.

```python
import inductiva

# Allocate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="n2-highcpu-2",
    spot=True
)

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST(
    version="4.0.2")

# Run simulation
task = openfast.run(
    input_dir="input_files",
    commands=[
        "openfast 5MW_OC4Semi_WSt_WavesWN/"
        "5MW_OC4Semi_WSt_WavesWN.fst"],
    on=cloud_machine)

task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()

```



As mentioned before, there is no point in trying to run this
simulation of a machine with many cores. So, we chose running
the sumulation on a `n2-highcpu-2` VM, that has 2 virtual CPUs, 
(that maps to 1 physical core under the hood). What is 
interesting is that `n2-highcpu-2` VM we chose is one of the least
expensive machines you can possibly get via Google Cloud: it costs
a mere 0.0081 US$/hour when running in spot mode.

So, since this is a relatively short simulation (30 seconds), running
it end-to-end on a n2-highcpu-2 costs only 0.0001 US$. 

## From 1 to 100
Inductiva does not help run one OpenFast simulation faster, but
it helps you run many simulations in parallel. So, let's assume
that you need to understand the impact of changing a certain 
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
at 1 meter steps. This means we are going to run 41 simulations.
But we are going to run these 41 simulation in parallel.

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


@Paulo: vou deixar isto agora para ti.
Ha aqui alguns detalhes:

1) convinha garantir que as tarefas passariam a ter uma
label mais explicita tipo "5MW_OC4Semi_WSt_WavesWN at WtrDpth = X".
Isto ja eh possivel? Se nao for, penso que deviamos fazer o push
por isso asap.

2) Mesmo que nao seja, convinha que houvesse forma de organizar 
as tasks em grupos / pastas / projectos. Isso ja eh possivel?
Se sim, vamos fazer isso ja nesse exemplo. Idealmente, seria implicito,
algo deste genero:

task = openfast.run(
    input_dir=target_dir,
    commands=["openfast 5MW_OC4Semi_WSt_WavesWN.fst"], 
    on=machine_group,
    task_set/project="Water Depth Exploration")

e nao explicito, em que temos mesmo de ter uma class Project, etc...

3) Nao sei muito bem como fazer um blocking para todas as tarefas
do grupo. Ou seja, esperar ate todas estarem completas. Imagino
que as coisas que foram feitas para os benchmarks ja permitam fazer
isso. Eu implementei manualmente assim:


import inductiva
import time


# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    machine_type="n2-highcpu-2",
    spot=True,
    num_machines=20
)
machine_group.start()

my_tasks = []
for depth in range(190, 210):

    print(f"Preparing files for depth = {depth}")
    target_dir = f"variations/params_for_depth_{depth}"

    inductiva.TemplateManager.render_dir(
        source_dir="openfast-5MW_OC4Semi_WSt_WavesWN_templated",
        target_dir=target_dir,
        overwrite=True,
        water_depth=depth)

    # Initialize OpenFAST simulator
    openfast = inductiva.simulators.OpenFAST()

    task = openfast.run(
        input_dir=target_dir,
        commands=["openfast 5MW_OC4Semi_WSt_WavesWN.fst"],
        on=machine_group)

    my_tasks.append(task)

while(not all([x.is_terminal() for x in my_tasks])):
    print("Waiting for ALL tasks to finish")
    print(f"Finished: {sum([x.is_terminal() for x in my_tasks])}")
    time.sleep(5)


machine_group.terminate()


Mas isto deveria ser muito facil.
Sera que podes ver como eh que isto seria feito de uma
forma elegante e ja integrada?
