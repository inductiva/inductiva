"""Calculix example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "calculix-input-example.zip",
    unzip=True)

# Initialize the Simulator
calculix = inductiva.simulators.Calculix( \
    version="2.22")

# Run simulation with config files in the input directory
task = calculix.run( \
    input_dir=input_dir,
    sim_config_filename="hueeber3.inp",
    on=cloud_machine,
    n_vcpus=4)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
