"""OpenFOAM ESI example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-esi-input-example.zip",
    unzip=True)

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM(distribution="esi")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir,
                    bash_script="./Allrun",
                    on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
