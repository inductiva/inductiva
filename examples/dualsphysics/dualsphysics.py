"""DualSPHysics example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("a3-highgpu-8g")

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir="path/to/my/DualSPHysics/files",
                        shell_script="my_config_file.sh",
                        on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
