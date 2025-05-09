"""OpenTelemac example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
telemac2d = inductiva.simulators.OpenTelemac( \
    version="8p4r0")

#  List of commands to run
commands = [
    # List the OpenTelemac commands
    # you wish to execute
]

# Run simulation
task = telemac2d.run(\
    input_dir="/Path/to/opentelemac-input-example",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
