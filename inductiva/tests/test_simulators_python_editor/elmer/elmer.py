"""Elmer example"""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import Elmer
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "elmer-input-example.zip",
    unzip=True)

# Initialize the Simulator
elmer = Elmer( \
    version="9.0")

# Run simulation
task = elmer.run( \
    input_dir=input_dir,
    commands=[
        "ElmerGrid 1 2 winkel.grd",
        "ElmerSolver case.sif -ipar 2 1 1"],
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
