"""FVCOM example"""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fvcom-input-example.zip",
    unzip=True)

# Initialize FVCOM Simulator
# Check available versions with the cli command "inductiva simulators list"
fvcom = inductiva.simulators.FVCOM(version="5.1.0")

# Run simulation with config files in the input directory
task = fvcom.run(input_dir=input_dir,
                 working_dir="run/",
                 case_name="tst",
                 n_vcpus=1,
                 debug=7,
                 on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
