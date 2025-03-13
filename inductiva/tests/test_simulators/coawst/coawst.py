"""COAWST Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "coawst-input-example.zip",
    unzip=True)

# Initialize the Simulator
coawst = inductiva.simulators.COAWST( \
    version="3.8")

# Run simulation
task = coawst.run( \
    input_dir=input_dir,
    sim_config_filename="ocean_ducknc.in",
    build_coawst_script="build_coawst.sh",
    n_vcpus=4,
    use_hwthread=True,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
