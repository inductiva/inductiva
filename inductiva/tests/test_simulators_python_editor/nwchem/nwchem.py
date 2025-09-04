"""NWChem example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import NWChem
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "nwchem-input-example.zip",
    unzip=True)

# Initialize the Simulator
nwchem = NWChem(version="7.2.3")

# Run simulation
task = nwchem.run( \
    input_dir=input_dir,
    sim_config_filename="h2o_sp_scf.nw",
    n_vcpus=1,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
