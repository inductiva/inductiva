"""Octopus example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "octopus-input-example.zip",
    unzip=True)

octopus = inductiva.simulators.Octopus( \
    version="16.1")

task = octopus.run( \
    input_dir=input_dir,
    commands=["octopus"],
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
