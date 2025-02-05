"""GX example."""
import inductiva

# Allocate machine
gpu_machine_group = inductiva.resources.MachineGroup("a3-highgpu-8g")

# Initialize the Simulator
gx = inductiva.simulators.GX()

# Run simulation with config files in the input directory
task = gx.run(input_dir="/Path/to/My/GX/Files",
              sim_config_filename="my_config_file.in",
              on=gpu_machine_group)

# Wait for the simulation to finish and download the results
task.wait()
gpu_machine_group.terminate()

task.download_outputs()
