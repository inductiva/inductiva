"""DualSPHysics example."""
import inductiva

# Allocate Google cloud machine
cloud_machine_gpu = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="a3-highgpu-1g")

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run( \
    input_dir="path/to/my/DualSPHysics/files",
    shell_script="my_config_file.sh",
    on=cloud_machine_gpu)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine_gpu.terminate()

task.download_outputs()
