""" OpenSees example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees()

# Run simulation with config files in the input directory
task = opensees.run(input_dir="/Path/to/My/OpenSees/Files",
                    sim_config_filename="my_config_file.tcl",
                    on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
