# Reducing Data Generation in Simulations

Running simulations often results in the creation of huge amounts of data. While
some of this data is necessary for latter analysis, other data might be temporary
or unnecessary, contributing to unnecessary storage and transfer costs.

Below are some steps that can help you reduce the amount of data generated by a
simulation and, ultimately, lower your storage and transfer expenses.

## Step 1: Configure the Simulator to Generate Less Data

Many simulators allow users to configure how much data is generated during a
simulation. By adjusting these settings, you can potentially reduce the amount
of unnecessary data that is created. Consider the following options:

- **Limit the outputs**: Check if the simulator allows you to specify the
outputs you need. For example, you may only need results at certain time steps
or specific variables. Limiting the outputs can greatly reduce the volume of
generated data.
- **Adjust output frequency**: If the simulator writes results at frequent
intervals, consider reducing the frequency of these writes.

By configuring your simulator to generate only the data you truly need, you can
significantly reduce the overall data footprint of your simulation.

## Step 2: Clean Up After the Simulation

Once your simulation is completed, you might end up with a large number of
temporary or intermediate files that are no longer needed. Running a cleanup
after the simulation can free up valuable storage space.

All integrated simulators allow you to pass an argument `on_finish_cleanup` that can be
a shell script or a list of commands that will run once your simulation is finished.
You should use this argument to delete any temporary or unwanted files to save
on storage costs.

### Example: OpenFOAM

OpenFOAM, creates a separate folder for each process running the simulation.
Each process generates data specific to that run, but in the end, OpenFOAM
combines all of the process data into a single directory. Therefore, the process
folders themselves may no longer be needed after the simulation completes.

You can pass, using the `on_finish_cleanup` argument, a shell script that deletes
all the temporary folders that are not needed once the simulation ends.

> Take a look into our recipes to see how this is done [here](https://inductiva.ai/guides/how-it-works/recipes/storage-related/index).

Keep this tip in mind when running your next simulation to reduce the amount of
data generated, keeping your [storage and transfer costs](https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost) in check.

```{banner_small}
:origin: how_it_works_save_storage
```