"""FDS example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fds-input-example.zip",
    unzip=True)

fds = inductiva.simulators.FDS()

task = fds.run(input_dir=input_dir,
               sim_config_filename="mccaffrey.fds",
               n_vcpus=1,
               on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
