"""OpenFOAM Foundation example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-8")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-input-example.zip",
    unzip=True)

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM( \
    distribution="foundation")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir,
                    shell_script="./Allrun",
                    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()
