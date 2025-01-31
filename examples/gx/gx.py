"""GX example."""
import inductiva

# Instantiate machine group
gpu_machine_group = inductiva.resources.MachineGroup("g2-standard-24")

gx = inductiva.simulators.GX()

task = gx.run(input_dir="/Path/to/My/GX/Files",
              sim_config_filename="my_config_file.in",
              on=gpu_machine_group)

task.wait()
gpu_machine_group.terminate()

task.download_outputs()

task.print_summary()
