"""FVCOM example"""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    machine_type="c3d-standard-180",
    provider="GCP")

# Initialize the Simulator
fvcom = inductiva.simulators.FVCOM()

# Run simulation with config files in the input directory
task = fvcom.run(input_dir="path/to/my/fvcom/files",
                 working_dir="run/",
                 case_name="my_case_name",
                 on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
