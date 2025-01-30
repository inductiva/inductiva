"""OpenFOAM ESI example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi", version="2412")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir="/path/to/my/openfoam/files",
                    shell_script="./Allrun",
                    on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
