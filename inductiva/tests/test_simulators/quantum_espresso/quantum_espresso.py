"""Quantum ESPRESSO example."""
import inductiva
from inductiva.commands import MPIConfig, Command

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
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

# Initialize QuantumEspresso simulator
qe = inductiva.simulators.QuantumEspresso( \
    version="7.4.1")

# Run simulation
task = qe.run( \
    input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
