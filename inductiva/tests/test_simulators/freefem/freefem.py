"""FreeFEM example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "freefem-input-example.zip",
    unzip=True)

# Initialize the Simulator
freefem = inductiva.simulators.FreeFEM( \
    version="4.15")

# Run simulation with config files in the input directory
task = freefem.run( \
    input_dir=input_dir,
    commands=[
        "FreeFem++ -nw mycode.edp"],
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
