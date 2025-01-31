"""Quantum ESPRESSO example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize QuantumEspresso simulator
qe = inductiva.simulators.QuantumEspresso()

my_qe_command = [
    # here you list the qe command you wish to execute
]

# Run simulation
task = qe.run(input_dir="/path/to/my/quantumEspresso/files",
              commands=my_qe_command,
              on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
