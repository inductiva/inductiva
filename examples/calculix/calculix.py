"""Calculix example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
calculix = inductiva.simulators.Calculix( \
    version="2.22")

# Run simulation with config files in the input directory
task = calculix.run( \
    input_dir="path/to/my/calculix/files",
    sim_config_filename="my_config_file.inp",
    on=cloud_machine,
    n_vcpus=4)

task.wait()
cloud_machine.terminate()

task.download_outputs()
