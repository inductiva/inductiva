"""CaNS example"""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cans-input-example.zip",
    unzip=True)

# Initialize the Simulator
cans = inductiva.simulators.CaNS(version="2.3.4")

# Run simulation with config files in the input directory
task = cans.run(input_dir=input_dir,
                sim_config_filename="input.nml",
                on=machine_group,
                n_vcpus=4)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
