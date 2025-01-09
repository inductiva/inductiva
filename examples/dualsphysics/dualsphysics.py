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

commands = [
    "gencase config flow_cylinder -save:all",
    "dualsphysics flow_cylinder flow_cylinder -dirdataout data -svres",
    ("partvtk -dirin flow_cylinder/data -savevtk flow_cylinder/PartFluid "
     "-onlytype:-all,+fluid")
]

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir=input_dir,
                        commands=commands,
                        on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
