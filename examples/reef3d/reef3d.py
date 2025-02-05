"""Reef3D example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize simulator
reef3d = inductiva.simulators.REEF3D()

# Run simulation
task = reef3d.run(input_dir="/path/to/my/reef3d/files", on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
