"""Delft3DFM example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "delft3dfm-input-example.zip",
    unzip=True)

delft3dfm = inductiva.simulators.Delft3DFM(version="2024.03")

# List of commands to run
commands = ["dflowfm --autostartstop f34.mdu"]

task = delft3dfm.run( \
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
