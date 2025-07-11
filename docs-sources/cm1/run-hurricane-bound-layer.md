# LES of Hurricane Bound Layer 🌪

In this tutorial, we demonstrate how to simulate a simplified hurricane boundary
layer (HBL) using Large Eddy Simulation (LES) on the Inductiva platform. This
setup follows the methodology of Bryan et al. (2017) **[1]**.

The hurricane boundary layer is the lowest portion of the atmosphere in a
tropical cyclone where turbulence, surface friction, and radial inflow are
dominant. This LES setup simulates the internal dynamics of the HBL without full
coupling to large-scale weather systems.

We will also demonstrate Inductiva’s ability to efficiently scale this use case,
starting with a cloud machine equivalent to a typical laptop and then scaling up
to more powerful machines.

## Prerequisites

The required files may be downloaded [here](https://storage.googleapis.com/inductiva-api-demo-files/cm1-les-hurr-bound-layer.zip).
Note that manual download is not necessary since the run script will donwload
it.

## Key Considerations

TODO

## Run the Simulation

Below is the code required to run the hurricane boundary layer use case with
the Inductiva API.

The simulation will then be sent to a `c2d-higcpu-16` virtual machine from
Google Cloud, equipped with 16 vCPU and 32 GB of RAM. This machine is
equivalent to a standard working laptop.

```python
"""LES of Hurricane Bound Layer."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True,
)

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cm1-les-hurr-bound-layer.zip",
    unzip=True,
)

# Initialize the Simulator
cm1 = inductiva.simulators.CM1(version="21.1")

# Run simulation with config files in the input directory
task = cm1.run(
    input_dir=input_dir,
    sim_config_filename="namelist.input",
    on=cloud_machine,
)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

Copy and paste it into a file named `run.py` in your working directory and
execute it by running:

````
python run.py
````

The script will take about **1 hour and 36 minutes** to run. In the end, you
should see something like:

```
Task status: Success

Timeline:
	Waiting for Input         at 09/07, 14:08:46      0.676 s
	In Queue                  at 09/07, 14:08:47      58.479 s
	Preparing to Compute      at 09/07, 14:09:45      4.285 s
	In Progress               at 09/07, 14:09:49      5772.233 s
		└> 5771.837 s      /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus cm1.exe namelist.input
	Finalizing                at 09/07, 15:46:02      6.8 s
	Success                   at 09/07, 15:46:08      

Data:
	Size of zipped output:    617.07 MB
	Size of unzipped output:  1.21 GB
	Number of output files:   776

Estimated computation cost (US$): 0.15 US$
```

## Scaling Up

| Machine Type     | Execution Time | Estimated Cost (USD) | Speedup   |
|------------------|----------------|----------------------|-----------|
| c2d-highcpu-16   | 1h 36min       | $0.14                | Reference |
| c4-highcpu-16    | 1h 6min        | $0.33                | 1.45×     |
| c2d-highcpu-112  | 17 min         | $0.17                | 5.65×     |
| c3d-highcpu-360  | 9 min          | $0.51                | 10.67×    |

## References

[1] [Simple hurricane boundary layer (HBL)](https://github.com/george-bryan/CM1/tree/333342b50c85577450868280c2d1cbeff90e2f89/run/config_files/les_HurrBoundLayer)
