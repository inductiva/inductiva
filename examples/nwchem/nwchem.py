"""NWChem example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
nwchem = inductiva.simulators.NWChem()

#Run simulation
task = nwchem.run(input_dir="/path/to/my/nwchem/files",
                  sim_config_filename="h2o_sp_scf.nw",
                  n_vcpus=90,
                  on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
