"""WAVEWATCH III example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "wavewatch3-input-example.zip",
    unzip=True)

# Initialize the Simulator
wavewatch3 = inductiva.simulators.WaveWatch3( \
    version="11-2024")

# Define the commands to run the simulation
commands = [ \
    "ww3_grid",
    "ww3_prep",
    "ww3_shel"]

# Run the simulation
task = wavewatch3.run( \
    input_dir=input_dir,
    custom_switch="switch_PR3_UQ_MPI",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
