"""AMR-Wind example"""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind()

# Run simulation with config files in the input directory
task = amr_wind.run( \
    input_dir="path/to/my/amr-wind/files",
    sim_config_filename="my_config_file.inp",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
