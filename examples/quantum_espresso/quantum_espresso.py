"""Quantum ESPRESSO example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize QuantumEspresso simulator
qe = inductiva.simulators.QuantumEspresso()

my_qe_command = [
    # List the QE commands you wish to execute
]

# Run simulation
task = qe.run(input_dir="/path/to/my/quantumEspresso/files",
              commands=my_qe_command,
              on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
