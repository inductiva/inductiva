""" SNL SWAN example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

# Initialize the Simulator
snl_swan = inductiva.simulators.SNLSWAN()

# Run simulation with config files in the input directory
task = snl_swan.run( \
    input_dir="/Path/to/My/Snl-Swan/Files",
    sim_config_filename="my_config_file.swn",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
