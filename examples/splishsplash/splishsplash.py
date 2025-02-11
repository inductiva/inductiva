"""SPlisHSPlasH example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    machine_type="c3d-standard-180",
    provider="GCP")

# Initialize the Simulator
splishsplash = inductiva.simulators.SplishSplash()

# RRun simulation with config files in the input directory
task = splishsplash.run(input_dir="/path/to/my/splishsplash/files",
                        sim_config_filename="my_config_file.json",
                        on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
