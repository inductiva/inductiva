"""OpenFOAM Foundation example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-8")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-input-example.zip",
    unzip=True)

# Run the Allrun script
commands = ["bash ./Allrun"]

# Or set the commands manually (Examples below)

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="foundation")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir, commands=commands, on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
