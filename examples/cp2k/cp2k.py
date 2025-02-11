""" CP2K example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    machine_type="c3d-standard-180",
    provider="GCP")

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K()

# Run simulation with config files in the input directory
task = cp2k.run(input_dir="/Path/to/My/cp2k/Files",
                sim_config_filename="my_config_file.inp",
                on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
