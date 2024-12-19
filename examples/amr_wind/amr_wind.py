"""AMR-Wind example"""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "amr-wind-input-example.zip",
    unzip=True)

# Initialize AmrWind Simulator
# Check available versions with the cli command "inductiva simulators list"
amr_wind = inductiva.simulators.AmrWind(version="1.4.0")

# Run simulation with config files in the input directory
task = amr_wind.run(input_dir=input_dir,
                    sim_config_filename="abl_amd_wenoz.inp",
                    on=machine_group,
                    n_vcpus=4)

task.wait()
task.download_outputs()

machine_group.terminate()
