"""OpenSees example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openseespy-input-example.zip",
    unzip=True)

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="python",
    version="3.7.1")

# Run simulation
task = opensees.run( \
    input_dir=input_dir,
    sim_config_filename="example_mpi_paralleltruss_explicit.py",
    use_hwthread=True,
    n_vcpus=4,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
