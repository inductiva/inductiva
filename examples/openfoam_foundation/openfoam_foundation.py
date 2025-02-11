"""OpenFOAM Foundation example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    machine_type="c3d-standard-180",
    provider="GCP")

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="foundation",
                                         version="12")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir="/path/to/my/openfoam/files",
                    shell_script="./Allrun",
                    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
