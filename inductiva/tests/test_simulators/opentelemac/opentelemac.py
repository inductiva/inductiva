"""Opentelemac example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-16")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "opentelemac-input-example.zip",
    unzip=True)

opentelemac = inductiva.simulators.OpenTelemac()

# List of commands to run
commands = ["telemac2d.py t2d_malpasset-fine.cas --ncsize=16"]

task = opentelemac.run( \
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
