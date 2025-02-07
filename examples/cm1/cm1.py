""" CM1 example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator with the desired mode ("mpi" or "openmp")
cm1 = inductiva.simulators.CM1(mode="mpi")

# Run simulation with config files in the input directory
task = cm1.run(input_dir="/Path/to/My/cm1/Files",
               sim_config_filename="my_config_file.inp",
               on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
