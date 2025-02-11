"""Gx example."""
import inductiva

# Instantiate machine group
gpu_cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="g2-standard-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "gromacs-input-example.zip",
    unzip=True)

gx = inductiva.simulators.GX()

task = gx.run( \
    input_dir=input_dir,
    sim_config_filename="itg_w7x_adiabatic_electrons.in",
    on=gpu_cloud_machine)

task.wait()
gpu_cloud_machine.terminate()

task.download_outputs()

task.print_summary()
