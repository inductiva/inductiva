""" CP2K example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K(
    version="2025.1")

# Run simulation with config files in the input directory
task = cp2k.run( \
    input_dir="/Path/to/My/cp2k/Files",
    sim_config_filename="my_config_file.inp",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
