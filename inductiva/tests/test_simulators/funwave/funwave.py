"""FUNWAVE example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Download input files and store them in a directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "funwave-input-example.zip",
    unzip=True)

# Initialize the Simulator
funwave = inductiva.simulators.FUNWAVE( \
    version="3.6")

# Run simulation
task = funwave.run( \
    input_dir=input_dir,
    sim_config_filename="input.txt",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
