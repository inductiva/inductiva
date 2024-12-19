"""OpenFOAM ESI example."""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-esi-input-example.zip",
    unzip=True)

# Initialize OpenFOAM Simulator
# Check available versions with the cli command "inductiva simulators list"
openfoam = inductiva.simulators.OpenFOAM(distribution="esi", version="2406")

# Run simulation with config files in the input directory
task = openfoam.run(input_dir=input_dir,
                    shell_script="./Allrun",
                    on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
