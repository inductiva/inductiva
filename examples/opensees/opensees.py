""" OpenSees example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees()

# Run simulation with config files in the input directory
task = opensees.run( \
    input_dir="/Path/to/My/OpenSees/Files",
    sim_config_filename="my_config_file.tcl",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
