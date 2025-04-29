"""COAWST Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
coawst = inductiva.simulators.COAWST( \
    version="3.8")

# Run simulation
task = coawst.run( \
    input_dir="/Path/to/My/COAWST/Files",
    sim_config_filename="my_config_file.in",
    build_coawst_script="build_coawst.sh",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()
