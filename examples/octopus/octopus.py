""" Octopus example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
octopus = inductiva.simulators.Octopus( \
    version="16.1")

# Run simulation with config files in the input directory
task = octopus.run( \
    input_dir="/Path/to/My/octopus/Files",
    commands=["octopus"],
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
