"""DualSPHysics example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir="path/to/my/DualSPHysics/files",
                        shell_script="my_config_file.sh",
                        on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
