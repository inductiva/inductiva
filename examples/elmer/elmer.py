"""Elmer example"""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
elmer = inductiva.simulators.Elmer()

my_elmer_command = [
    # List the elmer commands you wish to execute
]

# Run simulation with config files in the input directory
task = elmer.run( \
    input_dir="path/to/my/elmer/files",
    commands=my_elmer_command,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
