"""Reef3D example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

# Initialize simulator
reef3d = inductiva.simulators.REEF3D()

# Run simulation
task = reef3d.run(input_dir="/path/to/my/reef3d/files", on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
