"""OpenFOAM ESI example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi", version="2412")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir="/path/to/my/openfoam/files",
                    shell_script="./Allrun",
                    on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()