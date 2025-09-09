"""OpenSees example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import OpenSees
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "openseespy-input-example.zip",
    unzip=True)

# Initialize the Simulator
opensees = OpenSees( \
    interface="python",
    version="3.7.1")

# Run simulation
task = opensees.run( \
    input_dir=input_dir,
    sim_config_filename="example_mpi_paralleltruss_explicit.py",
    use_hwthread=True,
    n_vcpus=4,
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
