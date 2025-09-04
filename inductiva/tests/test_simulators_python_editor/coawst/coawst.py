"""COAWST Simulation."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import COAWST
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "coawst-input-example.zip",
    unzip=True)

# Initialize the Simulator
coawst = COAWST(version="3.8")

# Run simulation
task = coawst.run( \
    input_dir=input_dir,
    sim_config_filename="ocean_ducknc.in",
    build_coawst_script="build_coawst.sh",
    n_vcpus=4,
    use_hwthread=True,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
