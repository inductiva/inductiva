"""FVCOM example"""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")

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
                 on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
