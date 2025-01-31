"""GX example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

gx = inductiva.simulators.GX()

task = gx.run(input_dir="/Path/to/My/GX/Files",
              sim_config_filename="itg_w7x_adiabatic_electrons.in",
              on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
