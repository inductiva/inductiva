"""FreeFEM example"""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import FreeFEM
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "freefem-input-example.zip",
    unzip=True)

# Initialize the Simulator
freefem = FreeFEM( \
    version="4.15")

# Run simulation
task = freefem.run( \
    input_dir=input_dir,
    commands=[
        "FreeFem++ -nw mycode.edp"],
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
