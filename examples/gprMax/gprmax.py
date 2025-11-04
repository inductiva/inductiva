"""gprMax example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16")

# Initialize the Simulator
telemac2d = inductiva.simulators.gprMax( \
    version="v3.1.7", use_dev=True)

#  List of commands to run
commands = [
    "python -m gprmax input_file.in",
]

# Run simulation
task = telemac2d.run(\
    input_dir="/Path/to/gprmax-input-example",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
