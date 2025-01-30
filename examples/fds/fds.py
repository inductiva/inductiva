"""FDS example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

fds = inductiva.simulators.FDS()

task = fds.run(input_dir="path/to/my/fds/files",
               sim_config_filename="mccaffrey.fds",
               n_vcpus=90,
               on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
