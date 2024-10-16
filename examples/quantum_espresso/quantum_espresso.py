import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "qe-input-example.zip", unzip=True)

# List of commands to run
commands = [
    "pw.x -i Al_local_pseudo.in",
    "pw_openmp.x -i Al_qe_pseudo.in"
]

# Initialize QuantumEspresso simulator
qe = inductiva.simulators.QuantumEspresso()

# Run simulation 
task = qe.run(
        input_dir,
        commands=commands,
        n_vcpus=2,
        use_hwthread=False,
        on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()