"""AMR-Wind example"""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind()

# Run simulation with config files in the input directory
task = amr_wind.run(input_dir="path/to/my/amr-wind/files",
                    sim_config_filename="my_config_file.inp",
                    on=machine_group,
                    n_vcpus=90)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
