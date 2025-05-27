"""SFINCS example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
sfincs = inductiva.simulators.SFINCS( \
    version="2.2.1")

# Run simulation with configuration files in the input directory
task = sfincs.run( \
    input_dir="/path/to/my/sfincs/files",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
