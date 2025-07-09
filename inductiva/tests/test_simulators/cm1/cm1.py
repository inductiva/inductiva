"""CM1 example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cm1-input-example.zip",
    unzip=True)

# Initialize the Simulator
cm1 = inductiva.simulators.CM1( \
    version="21.1")

# Run simulation with config files in the input directory
task = cm1.run( \
    input_dir=input_dir,
    sim_config_filename="namelist.input",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
