"""DualSPHysics example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Download the configuration files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsphysics-input-example.zip",
    unzip=True)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir=input_dir,
                        shell_script="run.sh",
                        on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
