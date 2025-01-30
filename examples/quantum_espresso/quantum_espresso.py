"""Quantum ESPRESSO example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize QuantumEspresso simulator
qe = inductiva.simulators.QuantumEspresso()

# Run simulation
task = qe.run(input_dir="/path/to/my/quantumEspresso/files",
              commands=["pw_openmp.x -i Al_qe_pseudo.in"],
              on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
