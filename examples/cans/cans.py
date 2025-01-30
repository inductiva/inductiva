"""CaNS example"""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
cans = inductiva.simulators.CaNS()

# Run simulation with config files in the input directory
task = cans.run(input_dir="path/to/my/cans/files",
                sim_config_filename="input.nml",
                on=machine_group,
                n_vcpus=90)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
