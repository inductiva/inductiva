"""HEC-RAS example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c4-highcpu-4")

# Initialize the Simulator
hec_ras = inductiva.simulators.Hec( \
    distribution="ras")

# Specify the HEC-RAS commands you want to run, separated by commas
hec_ras_commands = [ \
    "RasGeomPreprocess Muncie.p04.tmp.hdf x04"]

# Run simulation
task = hec_ras.run( \
    input_dir="/path/to/my/hec-ras/files",
    commands=hec_ras_commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
