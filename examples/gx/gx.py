"""GX example."""
import inductiva

# Allocate Google cloud machine
cloud_machine_gpu = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="a3-highgpu-1g")

# Initialize the Simulator
gx = inductiva.simulators.GX()

# Run simulation with config files in the input directory
task = gx.run( \
    input_dir="/Path/to/My/GX/Files",
    sim_config_filename="my_config_file.in",
    on=cloud_machine_gpu)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine_gpu.terminate()

task.download_outputs()
