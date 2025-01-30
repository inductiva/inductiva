"""OpenFOAM ESI example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-esi-input-example.zip",
    unzip=True)

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir,
                    shell_script="./Allrun",
                    on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
