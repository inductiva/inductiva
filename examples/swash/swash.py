"""SWASH example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-input-example.zip",
    unzip=True)

# Initialize the Simulator
swash = inductiva.simulators.SWASH()
# or alternatively, to use a specific version of SWASH:
# swash = inductiva.simulators.SWASH(version="10.05")

# Run simulation with config files in the input directory
task = swash.run(input_dir=input_dir,
                 sim_config_filename="input.sws",
                 on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
