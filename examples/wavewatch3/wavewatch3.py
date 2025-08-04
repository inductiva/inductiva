"""WAVEWATCH III example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

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
    input_dir="/path/to/my/wavewatch3/files",
    switch="Ifremer2_pdlib",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
