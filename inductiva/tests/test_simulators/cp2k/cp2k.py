"""CP2K example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-highcpu-4")

# Download input files and store them in a directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cp2k-input-example.zip",
    unzip=True)

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K( \
    version="2025.1")

# Run simulation
task = cp2k.run( \
    input_dir=input_dir,
    sim_config_filename="Ac.inp",
    use_hwthread=True,
    n_vcpus=4,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
