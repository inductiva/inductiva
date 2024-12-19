"""Reef3D example."""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip",
    unzip=True)

# Initialize REEF3D Simulator
# Check available versions with the cli command "inductiva simulators list"
reef3d = inductiva.simulators.REEF3D(version="24.02")

task = reef3d.run(input_dir=input_dir, on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
