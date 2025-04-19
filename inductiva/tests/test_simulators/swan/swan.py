"""SWAN example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swan-input-example.zip", True)

# Initialize the Simulator
swan = inductiva.simulators.SWAN( \
    version="41.45")

# Run simulation with config files in the input directory
# Uses swanrun by default
task = swan.run( \
    input_dir=input_dir,
    sim_config_filename="a11refr.swn",
    on=cloud_machine)

# Uses swanrun
# task = swan.run(input_dir=input_dir,
#                 sim_config_filename="a11refr.swn",
#                 command="swanrun",
#                 on=cloud_machine)

# Uses swan.exe. Note: the simulation file must be called INPUT
# task = swan.run(input_dir=input_dir,
#                 command="swan.exe",
#                 on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
