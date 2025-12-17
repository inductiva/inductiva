"""Delft3DFM example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
delft3dfm = inductiva.simulators.Delft3DFM( \
    version="2024.03")

# List of commands to run
commands = [
    "dflowfm --autostartstop f34.mdu",
]

# Run simulation
task = delft3dfm.run( \
    input_dir="path/to/my/delft3dfm/files",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
