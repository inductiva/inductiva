"""Quantum ESPRESSO example."""
import inductiva
from inductiva.commands import MPIConfig, Command

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "qe-input-example.zip",
    unzip=True)

mpi_config = MPIConfig(version="4.1.6", np=2, use_hwthread_cpus=False)

# List of commands to run
commands = [
    Command("pw.x -i Al_local_pseudo.in", mpi_config=mpi_config),
    # openMP command should not be used with MPI
    "pw_openmp.x -i Al_qe_pseudo.in"
]

# Initialize QuantumEspresso simulator
qe = inductiva.simulators.QuantumEspresso()

# Run simulation
task = qe.run(input_dir, commands=commands, on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
