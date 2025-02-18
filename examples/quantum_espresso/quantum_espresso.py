"""Quantum ESPRESSO example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

# Initialize QuantumEspresso simulator
qe = inductiva.simulators.QuantumEspresso()

my_qe_command = [
    # List the QE commands you wish to execute
]

# Run simulation
task = qe.run( \
    input_dir="/path/to/my/quantumEspresso/files",
    commands=my_qe_command,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
