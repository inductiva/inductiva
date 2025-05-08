"""OpenTelemac example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(provider="GCP",
                                                 machine_type="c2d-highcpu-16")

# Initialize the Simulator
telemac2d = inductiva.simulators.OpenTelemac( \
    version="8p4r0")

#  List of commands to run
commands = [
    # Run the simulation using 32 cores
    "telemac2d.py t2d_malpasset-fine.cas --ncsize=16",
    # Convert the results to VTK format
    "converter.py srf2vtk r2d_malpasset-fine.slf t2d_malpasset.vtk",
]

# Run simulation
task = telemac2d.run(input_dir="/Path/to/opentelemac-input-example",
                     commands=commands,
                     on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
