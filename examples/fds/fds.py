import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup('c2-standard-4')
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fds-input-example.zip", unzip=True)

fds = inductiva.simulators.FDS()

task = fds.run(input_dir=input_dir,
               sim_config_filename="mccaffrey.fds",
               post_processing_filename="mccaffrey.ssf",
               n_vcpus=1,
               on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()