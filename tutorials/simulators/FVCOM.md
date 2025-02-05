In this guide, we will walk you through setting up and running FVCOM, one 
of the built-in simulators available through the Inductiva API. 

We will cover:

- An overview of how we compiled FVCOM (fvcom and fvcom_estuary binaries.)
- Generating a valid namelist file for configuring your simulations.
- Example code to help you get started with simulations.

# FVCOM

FVCOM (Finite Volume Community Ocean Model)
is a 3D hydrodynamic model specifically designed for simulating coastal 
and ocean dynamics. It uses an unstructured grid and finite-volume methods, 
making it highly adaptable for modeling complex coastlines, estuaries, 
and bathymetry.

FVCOM excels at simulating ocean circulation, tides, and coastal processes, 
providing high-resolution outputs for water currents, temperature, salinity, 
and ecosystem interactions. It’s a valuable tool for studying marine 
ecosystems, coastal environments, and the impacts of climate change on 
ocean systems.

## Running FVCOM Simulations

We have compiled two versions:

- fvcom: The standard version for general use.
- fvcom_estuary: Configured to run the Estuary test case included with FVCOM.

Users can check the compilation modules used in each binary in the 
respective `make.inc` files located at `/make.inc` and `/make_estuary.inc`.

If you encounter issues with the input namelist file, the following Python 
script will generate a valid `.nml` file in your working directory:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fvcom-input-example.zip", unzip=True)

# Initialize the Simulator
fvcom = inductiva.simulators.FVCOM()

# Run simulation with config files in the input directory
task = fvcom.run( input_dir=input_dir,
                  working_dir="run/",
                  create_namelist="tst",
                  n_vcpus=1,
                  on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()
```

## Example Code

Below is an example where we run a simple FVCOM test scenario to verify 
that the simulator is working correctly, using **1 MPI process** and a **debug level of 7**:

```{literalinclude} ../../examples/fvcom/fvcom.py
:language: python
```

**Closing Notes**: If you encounter issues with the timezone argument in the `nml` file, please set it to `None` or `UTC` as a workaround. For more information, 
check the [bug report](https://github.com/FVCOM-GitHub/FVCOM/issues/27).

All simulation parameters were changed to make sure it runs in a reasonable time.
