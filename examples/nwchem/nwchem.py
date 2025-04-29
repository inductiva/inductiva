"""NWChem example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
nwchem = inductiva.simulators.NWChem()

# Run simulation with config files in the input directory
task = nwchem.run( \
    input_dir="/path/to/my/nwchem/files",
    sim_config_filename="my_config_file.nw",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
