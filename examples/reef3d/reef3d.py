"""Reef3D example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize simulator
reef3d = inductiva.simulators.REEF3D()

# Run simulation
task = reef3d.run(input_dir="/path/to/my/reef3d/files", on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
