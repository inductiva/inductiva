"""SCHISM example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

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
task.download_outputs()

machine_group.terminate()
