"""Quantum ESPRESSO example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import QuantumEspresso
from inductiva.utils import download_from_url
from inductiva.commands import MPIConfig, Command

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "qe-input-example.zip",
    unzip=True)

mpi_config = MPIConfig( \
    version="4.1.6",
    np=2,
    use_hwthread_cpus=False)

# List of commands to run
commands = [
    Command("pw.x -i Al_local_pseudo.in", mpi_config=mpi_config),
    # openMP command should not be used with MPI
    "pw_openmp.x -i Al_qe_pseudo.in"
]

# Initialize simulator
qe = QuantumEspresso(version="7.4.1")

# Run simulation
task = qe.run( \
    input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
