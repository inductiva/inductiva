"""HEC-RAS example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "hec-ras-input-example.zip",
    unzip=True)

# Initialize the Simulator
hec_ras = inductiva.simulators.Hec( \
    distribution="ras")

# Specify the HEC-RAS commands you want to run, separated by commas
hec_ras_commands = [
    "RasGeomPreprocess Muncie.p04.tmp.hdf x04",
    "mv Muncie.p04.tmp.hdf Muncie.p04.hdf",
    "python3 remove_HDF5_Results_Sed.py Muncie.p04.hdf",
    "RasUnsteady Muncie.p04.tmp.hdf x04",
    "mv Muncie.p04.tmp.hdf Muncie.p04.hdf", "RasSteady Muncie.r04"
]

# Run simulation
task = hec_ras.run( \
    input_dir=input_dir,
    commands=hec_ras_commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
