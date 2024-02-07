# Templating engine

The Inductiva API is all about enabling you to simulate at scale. As we have shown, with a few lines of Python code, you can send your simulations to MPI Clusters assembled from last-generation cloud hardware, letting you run much larger simulations than you would be able to using your local resources. Or you can spin up a large Machine Group, with dozens or hundreds of machines, and send a large number of simulations to be run on those machines in parallel. And such massive parallelism is precisely what you need when you are developing projects that require simulating a large number of variations of a certain base scenario. 

For example, suppose that for protecting a certain area on the coast, you are trying to find the best location and orientation for building a simple seawall, which can still have a number of possible variations. Suppose that you have 10 possible variations on the shape of the seawall, and you are considering 20 possible locations where you could build the wall. Also, you can direct the wall in 10 possible orientations. This design space includes a total of 2000 variations. For each configuration, you want to simulate the effectiveness of the wall under 25 different sea and wind conditions, which represent the extreme cases the wall is supposed to protect against. To completely explore the space of design and test it under all selected conditions, we would need to run 50000 simulations. 

This is quite a large number, and this is why Python scripting is great for this kind of challenge: you can build a set of “for loops'' that transverse all possible values for each parameter (shape, location, orientation, sea and wind conditions) and issue the corresponding simulation. Since the Inductiva API allows you to build large Machine Groups (e.g. 1000 machines) then you could potentially run 1000 of these variations in parallel and execute the 50000 simulation approximately in the time it would take to run 50 simulations in one machine (which we are going to assume it is a reasonable amount of time).

But, if each simulation is configured using a set of files, how do we programmatically change those simulation configuration files so that we can run 50000 *different* simulations, each one being a slight variation of another. 

This is where Inductiva’s templating mechanism comes into play. Templating allows you to start with a specific simulation file – your “base case” – containing fixed values for the parameters you wish to explore, and transform those fixed values into variables that you can now change programmatically from your Python code before you submit the simulation for remote execution. Let's illustrate the power of templating in a simple simulation case, from which you will be able to generalize to your own cases.

A simple example: Experimenting fluids with different properties

Suppose you want to study how fluids with different properties fall inside a container with a cubic shape. More specifically, you want to study the effect of fluid density and kinematic viscosity in the splashing on the walls of the cubic container. For that, you wish to experiment with 4 values of fluid density, and 4 values of kinematic viscosity. Overall, you will need to run 16 simulations. 

For the simulation of this scenation, you choose the SplishSplash simulator, and you start by coding your base case, where you assume the fluid is water. The corresponding simulation configuration file is the following:

@ivan 

Using the API, one would run the simulation on a XXX VM with this code:

@ivan 

This simulation takes about XX to run and the result looks like this (visualization produced with XXX):


Observe that the section of the simulation script that defines the properties of fluid is this one:

    "Materials": [
        {
            "id": "Fluid",
            "density0": XX,
            "viscosity": YY,
            "viscosityMethod": 6
        }



The density and kinematic viscosity shown above correspond to water. Now, we will edit this section and substitute each of the numeric values by descriptive label inside double curly brackets:

    "Materials": [
        {
            "id": "Fluid",
            "density0": {{ density }},
            "viscosity": {{ kinematic_viscosity }},
            "viscosityMethod": 6
        }


We can save this new file as ZZZ, and that’s it: we have a template. The labels we just assigned will now be treated as names of variables whose values we can programmatically assign using Python. 

We can use this template to run the simulation of fluids with different properties that are set programmatically. This is how we would do it for the fluid honey:

@ivan pls show just one simulation, not the 16 yet. Choose any other fluid you want.


The last step for generating the results we need for the study is to issue the 16 simulations from a single script by iterating over the set of admissible values for both parameters. Additionally, we can save time by issuing these 16 simulations in parallel. For that, we first need to create a machine group with 16 VMs:



Overall, the exploration of these 16 possibilities took xxx, approximately the same time needed to run only one simulation. And this was all done with minor changes to the simulation script that described the base case.

