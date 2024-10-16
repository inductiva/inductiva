import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip",
    unzip=True)

reef3d = inductiva.simulators.REEF3D()

task = reef3d.run(input_dir=input_dir, on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
