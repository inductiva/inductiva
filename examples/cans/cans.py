"""CaNS example"""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

# Initialize the Simulator
cans = inductiva.simulators.CaNS()

# Run simulation with config files in the input directory
task = cans.run(input_dir="path/to/my/cans/files",
                sim_config_filename="my_config_file.nml",
                on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
