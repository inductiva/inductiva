"""AMR-Wind example"""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind()

# Run simulation with config files in the input directory
task = amr_wind.run(input_dir="path/to/my/amr-wind/files",
                    sim_config_filename="my_config_file.inp",
                    on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()