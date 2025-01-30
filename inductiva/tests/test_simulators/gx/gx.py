"""Gx example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("g2-standard-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "gromacs-input-example.zip",
    unzip=True)

gx = inductiva.simulators.Gx()

task = gx.run(input_dir=input_dir,
              sim_config_filename="itg_w7x_adiabatic_electrons.in",
              on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
