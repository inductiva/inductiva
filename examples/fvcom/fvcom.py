"""FVCOM example"""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
fvcom = inductiva.simulators.FVCOM()

# Run simulation with config files in the input directory
task = fvcom.run(input_dir="path/to/my/fvcom/files",
                 working_dir="run/",
                 case_name="my_case_name",
                 n_vcpus=90,
                 debug=7,
                 on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
