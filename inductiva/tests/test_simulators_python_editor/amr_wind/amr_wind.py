"""AMR-Wind example"""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import AmrWind
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "amr-wind-input-example.zip",
    unzip=True)

# Initialize the Simulator
amr_wind = AmrWind( \
    version="3.4.1")

# Run simulation
task = amr_wind.run( \
    input_dir=input_dir,
    sim_config_filename="abl_amd_wenoz.inp",
    on=machine,
    n_vcpus=4)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
