"""CP2K example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import CP2K
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "cp2k-input-example.zip",
    unzip=True)

# Initialize the Simulator
cp2k = CP2K( \
    version="2025.1")

# Run simulation
task = cp2k.run( \
    input_dir=input_dir,
    sim_config_filename="Ac.inp",
    use_hwthread=True,
    n_vcpus=4,
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
