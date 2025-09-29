"""Elmer example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "elmer-input-example.zip",
    unzip=True)

# Initialize the Simulator
elmer = inductiva.simulators.Elmer( \
    version="9.0")

# Run simulation with config files in the input directory
task = elmer.run( \
    input_dir=input_dir,
    commands=[
        "ElmerGrid 1 2 winkel.grd",
        "ElmerSolver case.sif -ipar 2 1 1"],
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
