"""DualSPHysics example."""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

# Download the configuration files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsphysics-input-example.zip",
    unzip=True)

commands = [
    "gencase config flow_cylinder -save:all",
    "dualsphysics flow_cylinder flow_cylinder -dirdataout data -svres",
    ("partvtk -dirin flow_cylinder/data -savevtk flow_cylinder/PartFluid "
     "-onlytype:-all,+fluid")
]

# Initialize DualSPHysics Simulator
# Check available versions with the cli command "inductiva simulators list"
dualsphysics = inductiva.simulators.DualSPHysics(version="5.2.1")

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir=input_dir,
                        commands=commands,
                        on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
