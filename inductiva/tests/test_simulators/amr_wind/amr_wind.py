"""AMR-Wind example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "amr-wind-input-example.zip",
    unzip=True)

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind( \
    version="3.4.1")

# Run simulation with config files in the input directory
task = amr_wind.run( \
    input_dir=input_dir,
    sim_config_filename="abl_amd_wenoz.inp",
    on=cloud_machine,
    n_vcpus=4)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
