"""Gx example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
gx = inductiva.simulators.Gx()

task = gx.run(input_dir="/path/to/my/gx/files",
              sim_config_filename="itg_w7x_adiabatic_electrons.in",
              on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
