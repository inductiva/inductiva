"""SWMM example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-16")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swmm-input-example.zip",
    unzip=True)

# Initialize the Simulator
swmm = inductiva.simulators.SWMM(version="5.2.4")

#  List of commands to run
commands = [
    "runswmm model.inp model.rpt",
]

# Run simulation
task = swmm.run(\
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
