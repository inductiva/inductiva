"""SWASH example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-input-example.zip",
    unzip=True)

# Initialize the Simulator
swash = inductiva.simulators.SWASH(version="10.05")

# Run simulation with config files in the input directory
task = swash.run( \
    input_dir=input_dir,
    sim_config_filename="input.sws",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
