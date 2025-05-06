"""Delft3D example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
delft3d = inductiva.simulators.Delft3D( \
    version="6.04.00")

# List of commands to run
commands = ["mpirun -np 4 d_hydro.exe config_d_hydro.xml"]

# Run simulation
task = delft3d.run( \
    input_dir="path/to/my/delft3d/files",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
