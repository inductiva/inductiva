"""Gx example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("g2-standard-4")

gromacs = inductiva.simulators.Gx()

task = gromacs.run(input_dir="/Path/to/My/Gx/Files",
                   sim_config_filename="itg_w7x_adiabatic_electrons.in",
                   on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
