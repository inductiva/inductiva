"""OpenFAST example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST()

# Run simulation
task = openfast.run(input_dir="/path/to/my/openfast/files",
                    commands=["openfast IEA-15-240-RWT-Monopile.fst"],
                    on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
