"""XBeach example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach(version="1.24")

# Run simulation with configuration files in the input directory
task = xbeach.run( \
    input_dir="/path/to/my/xbeach/files",
    sim_config_filename="my_config_file.txt",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
