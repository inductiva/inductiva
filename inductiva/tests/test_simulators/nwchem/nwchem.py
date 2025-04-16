"""NWChem example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "nwchem-input-example.zip",
    unzip=True)

nwchem = inductiva.simulators.NWChem()

task = nwchem.run( \
    input_dir=input_dir,
    sim_config_filename="h2o_sp_scf.nw",
    n_vcpus=1,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
