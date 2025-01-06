"""SWAN example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swan-input-example.zip", True)

# Initialize the Simulator
swan = inductiva.simulators.SWAN()

# Run simulation with config files in the input directory
# Uses swanrun by default
task = swan.run(input_dir=input_dir,
                sim_config_filename="a11refr.swn",
                on=machine_group)

# Uses swanrun
# task = swan.run(input_dir=input_dir,
#                 sim_config_filename="a11refr.swn",
#                 command="swanrun",
#                 on=machine_group)

# Uses swan.exe. Note: the simulation file must be called INPUT
# task = swan.run(input_dir=input_dir,
#                 command="swan.exe",
#                 on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
task.download_outputs()

machine_group.terminate()

task.print_summary()
