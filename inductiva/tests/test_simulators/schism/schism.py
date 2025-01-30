"""SCHISM example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "schism-input-example.zip", True)

# Initialize the Simulator
schism = inductiva.simulators.SCHISM()

# Run simulation with config files in the input directory
task = schism.run(input_dir=input_dir,
                  n_vcpus=3,
                  num_scribes=2,
                  on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
