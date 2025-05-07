# Private docker usage


Inductiva provides a list of open-source simulation containers publicly available on [DockerHub - inductiva/kutu](https://hub.docker.com/r/inductiva/kutu). However, we understand that these may not cover all your simulation need, specially if you have application-specific requirements that involve running custom code as a pre-processing or post-processing stage.

To give users full control over the simulation environment, Inductiva supports Custom Images — allowing users to run simulations using their own Docker containers. The best part? These images don't need to be made public — Inductiva securely manages the storage and execution of your private containers.

## How does it work?

At Inductiva, we use Apptainer to run containerized simulations. Apptainer (formerly Singularity) is an open-source container platform tailored for high-performance computing (HPC).

Apptainer uses a specialized format called Singularity Image Format (SIF) — a single-file container image that's portable, secure, and easy to distribute.

When running simulations, Inductiva automatically converts Docker images into SIF format before execution. This enables the integration between user-defined Docker environments and our HPC infrastructure.

## How to create a Singularity Image Format (SIF) File?

To use your Docker image with Inductiva, you'll need to convert it to the SIF format. Fortunately, this is straightforward using the Inductiva CLI:

```bash
inductiva containers convert <DOCKER IMAGE NAME/ID/URL> <OUTPUT FILE>
```

- This command looks for the specified Docker image in your local Docker registry.
- It spins up a conversion container with Apptainer pre-installed.
- The container then converts the Docker image into a .sif file.

You can also reference remote DockerHub images by prefixing them with docker://, for example:

```bash
inductiva containers convert docker://username/my-private-image:latest private-image.sif
```

## Saving the SIF format into Inductiva Storage

To use the container in your simulations, the SIF file must be uploaded to your Inductiva storage.

While you can upload the .sif file manually, the easiest method is to convert and upload in a single step:

```bash
inductiva containers upload <DOCKER IMAGE NAME/ID/URL>
```

This command performs the following:

- Converts the Docker image to .sif format.
- Uploads the result directly to your private storage on Inductiva.

Once uploaded, the custom container is ready to be used via the Inductiva API for your simulations.


## Running Simulations with a Custom Container


Once your custom Docker image has been uploaded (in `.sif` format) to Inductiva storage, you can reference it directly in your simulation script using the `inductiva.simulators.CustomImage` class.

```python
import inductiva

# Retrieve your machine group (where the simulation will run)
mg = inductiva.resources.MachineGroup(...)

# Reference your uploaded SIF image in Inductiva Storage
inductiva_image = "inductiva://my-containers/private_image.sif"

# Create a simulator using your custom container image
simulator = inductiva.simulators.CustomImage(container_image=inductiva_image)

# Run the simulation
task = simulator.run(
    input_dir="./input_dir",                    
    commands=["echo Hello World!"],             
    on=mg                                       
)
```

Notes:

- The inductiva:// prefix is used to reference files stored in your Inductiva private storage.
- The commands list defines what will be executed inside the container once it runs.
- You can prepare your input files inside ./input_dir, which will be mounted inside the container.

Once the simulation is submitted, Inductiva will handle the environment setup, container execution, and output management automatically.