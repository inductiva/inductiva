"""FVCOM example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fvcom-input-example.zip",
    unzip=True)

# Initialize the Simulator
fvcom = inductiva.simulators.FVCOM()

# Run simulation with config files in the input directory
task = fvcom.run(input_dir=input_dir,
                 working_dir="run/",
                 case_name="tst",
                 n_vcpus=1,
                 debug=7,
                 on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
