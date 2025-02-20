# Advanced OpenFAST Tutorial

In this tutorial, we are going to show you how you can
use the Inductiva API to accelerate your OpenFast projects.

## The peculiar case of OpenFast

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
could you possibly have if you decide to send your OpenFast
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

<h3>🌊 OpenFAST Simulation Visualization</h3>
<p align="center">
  <img src="../_static/openfast_animation_30_fps.gif" alt="OpenFAST simulation visualization" width="700">
</p>

## Our use case: 5MW_OC4Semi_WSt_WavesWN

In this tutorial, we will show you how to do this using the 
"5MW_OC4Semi_WSt_WavesWN" example, available from the OpenFast 
[GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.
 
This example is an extension of the reference case described in 
["Definition of a 5-MW Reference Wind Turbine for Offshore
System Development"](https://www.nrel.gov/docs/fy09osti/38060.pdf).

All files needed are available from OpenFast GitHub so
let's get started.



<p align="center">
  <img src="./openfastAdvanced_static/single_turbine.jpg" alt="Openfast Simulation" width="300">
</p>

### Requirements: Setting up your files

Before we start running the simulation, we need to ensure that all required
files are properly set up. In this section, we'll go through the steps to
download, and prepare the necessary input files.

#### Step 1: Downloading you simulation files

We are going to run the `5MW_OC4Semi_WSt_WavesWN` case as
it is originally defined in the [GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.

We start by downloading the folders `5MW_OC4Semi_WSt_WavesWN` and the `5MW_Baseline` to our local input directory. This folders are located [here](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast).

After downloading the folders and moving them your `input_files` should look something like this:

- input_files
    - 5MW_Baseline
    - 5MW_OC4Semi_WSt_WavesWN

#### Step 2: Building `DISCON_OC3Hywind.dll`

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
Don't forget to paste this file in the 5MW_Baseline/ServoData folder.

You now have all the necessary files to run your simulation.

## Running our simulation once

For you to run this simulation you can execute the following python code.

```python
import inductiva

# Allocate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="n2-highcpu-2",
    spot=True)

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
interesting is that this `n2-highcpu-2` VM we chose is one of the least
expensive machines you can possibly get via Google Cloud: it costs
a mere 0.0081 US$/hour when running in spot mode.

So, since this is a relatively short simulation (29.968 seconds), running
it end-to-end on a n2-highcpu-2 costs only 0.00011 US$.

In the next part of this tutorial, we’ll take things to the next level by
running dozens of OpenFAST simulations in parallel using Inductiva,
demonstrating the true power of cloud-based scalability. Stay tuned!
