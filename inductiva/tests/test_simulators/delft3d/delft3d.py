"""Delft3D example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "delft3d-input-example.zip",
    unzip=True)

delft3d = inductiva.simulators.Delft3D(version="6.04.00")

# List of commands to run
commands = ["mpirun -np 4 d_hydro.exe config_d_hydro.xml"]

task = delft3d.run( \
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
