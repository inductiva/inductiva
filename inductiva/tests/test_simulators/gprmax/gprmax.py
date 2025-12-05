"""GprMax example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "gprmax-input-example.zip",
    unzip=True)

# Initialize the Simulator
gprmax = inductiva.simulators.GprMax(version="3.1.7")

#  List of commands to run
commands = [
    "python -m gprMax antenna_like_MALA_1200_fs.in",
]

# Run simulation
task = gprmax.run(\
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
