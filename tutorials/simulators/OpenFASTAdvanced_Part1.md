# OpenFAST Tutorial (Advanced)

In this tutorial, we are going to show you how you can
use the Inductiva API to accelerate your OpenFAST projects, by showing how to
run dozens of simulations in parallel.

<p align="center"><img src="../_static/openfast_animation_30_fps.gif" alt="OpenFAST simulation visualization" width="700"></p>

## The curious case of OpenFAST

OpenFAST is a pretty peculiar piece of software because
it runs only on a single thread. That means that OpenFAST
does not take full advantage of the extra compute capacity
offered by modern CPUs, which support dozens, or sometimes 
hundreds, of computational threads in parallel. In fact,
the most decisive factor impacting the performance of 
OpenFAST simulations is CPU clock frequency. 

Now, what is interesting about that is that, tipically, your 
500\$ home desktop has a much higher clock frequency than the 
20000\$ monster machines available in cloud centers. Your local
desktop has CPU with clock frequencies between 4 and 5GHz (but
can be overclocked to more than 5.5GHz), while cloud machines
have CPUs that will be running at around 3GHz (or less!).

But why is that?

Have you noticed how fast the ventilation 
fans of your desktop PC spin when you are running some long 
simulation? Now imagine how many fans and sophisticated
cooling systems you would need if you packed 100 CPUs like
yours in a volume the size of your kitchen fridge. Hardware
density in datacenter cabinets is so high that if all the
CPUs were running at speeds of 5GHz or more, the datacenter
would be impossible to cool down.

So, if that is the case, why would you even bother using
Inductiva to run your OpenFAST simulation? What advantage
could you possibly have if you decide to send your OpenFAST
simulation to a cloud machine using Inductiva? 

And the answer is: **hundreds of advantages!**

If you are running only 1 OpenFAST simulation, 
then you should run it on your desktop machine: it will be faster
there due to much higher clock speeds.  But if you need to run
hundreds or thousands of OpenFAST simulations, then you can use
Inductiva to spin up hundreds of very cheap cloud machines to 
run all those simulations in parallel instead of having to run
them sequentially in your machine! This will be hundreds of times
faster. And it is super easy (and cost-effective) with Inductiva.

## Our use case: 5MW_OC4Semi_WSt_WavesWN

In this tutorial, we will show you how to do this using the 
"5MW_OC4Semi_WSt_WavesWN" example, available on the OpenFAST 
[GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.
 
This example is an extension of the reference case described in 
["Definition of a 5-MW Reference Wind Turbine for Offshore
System Development"](https://www.nrel.gov/docs/fy09osti/38060.pdf).

All the necessary files are available from the OpenFAST GitHub, so let's get
started.

At the top of this page you will find a visualization of the simulation results.

In the following pages you will see the following topics:
- [Preparation](OpenFASTAdvanced_Part2.md):
    Download and prepare all the files needed to run our simulations.
- [Running the Simulation using Inductiva's API](OpenFASTAdvanced_Part3.md):
    Run two simulations and compare the results.
- [Running 50 Simulations in parallel - Templating](OpenFASTAdvanced_Part4.md):
    Modify our input files to support our templating system, for easy
    parameterization of some simulation variables.
- [Running 50 Simulations](OpenFASTAdvanced_Part5.md):
    Run 50 simulations in parallel.
- [Downloading The Results](OpenFASTAdvanced_Part6.md):
    Go through the simulation statuses and we also show you how to download the
    files of all 50 simulations to a specific folder for later post-processing
    on your local machine.

[Preparation](OpenFASTAdvanced_Part2.md)