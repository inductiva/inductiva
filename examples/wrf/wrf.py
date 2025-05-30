"""WRF Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
wrf = inductiva.simulators.WRF( \
    version="4.6.1")

# Run simulation
task = wrf.run( \
    input_dir="/Path/to/my/wrf/files",
    case_name="your_case_name",  # e.g., "em_real"
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()
