"""gprMax example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
gprmax = inductiva.simulators.GprMax(version="3.1.7")

#  List of commands to run
commands = [
    "python -m gprMax input_file.in",
]

# Run simulation
task = gprmax.run(\
    input_dir="/Path/to/gprmax-input-example",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
