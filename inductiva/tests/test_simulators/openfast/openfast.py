"""OpenFAST example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfastv4.0.2-input-example.zip",
    unzip=True)

# List of commands to run
commands = ["openfast 5MW_OC4Semi_WSt_WavesWN/"
            "5MW_OC4Semi_WSt_WavesWN.fst"]

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST( \
    version="4.0.2")

# Run simulation
task = openfast.run( \
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
