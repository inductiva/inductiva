"""SWMM example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c4-highcpu-2")

# Initialize the Simulator
swmm = inductiva.simulators.SWMM(version="5.2.4")

#  List of commands to run
commands = [
    "runswmm model.inp model.rpt",
]

# Run simulation
task = swmm.run(\
    input_dir="/Path/to/swmm-input-example",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
