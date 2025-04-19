"""SNLSWAN example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "snlswan-input-example.zip", True)

# Initialize the Simulator
snlswan = inductiva.simulators.SNLSWAN( \
    version="2.2")

# Run simulation with config files in the input directory
# Uses swanrun by default
task = snlswan.run( \
    input_dir=input_dir,
    sim_config_filename="input.swn",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
