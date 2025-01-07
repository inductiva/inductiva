"""CaNS example"""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cans-input-example.zip",
    unzip=True)

# Initialize the Simulator
cans = inductiva.simulators.CaNS()

# Run simulation with config files in the input directory
task = cans.run(input_dir=input_dir,
                sim_config_filename="input.nml",
                on=machine_group,
                n_vcpus=4)

task.wait()
task.download_outputs()

machine_group.terminate()

task.print_summary()
