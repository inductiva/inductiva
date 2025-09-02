"""SPlisHSPlasH example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip",
    unzip=True)

# Set simulation input directory
splishsplash = inductiva.simulators.SplishSplash( \
    version="2.13.0")

task = splishsplash.run( \
    input_dir=input_dir,
    sim_config_filename="config.json",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
