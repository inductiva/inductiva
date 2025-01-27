"""NWChem example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "nwchem-input-example.zip",
    unzip=True)

nwchem = inductiva.simulators.NWChem()

task = nwchem.run(input_dir=input_dir,
                  sim_config_filename="h2o_sp_scf.nw",
                  n_vcpus=1,
                  on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
