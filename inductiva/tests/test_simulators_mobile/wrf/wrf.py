"""WRF Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "wrf-input-example.zip",
    unzip=True)

# Initialize the Simulator
wrf = inductiva.simulators.WRF( \
    version="4.6.1")

# Run simulation
task = wrf.run( \
    input_dir=input_dir,
    init_commands=["./ideal.exe"],
    case_name="em_fire",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
