# How It Works?
Inductiva provides a collection of open-source simulation containers, 
publicly available through our Docker Hub repository, [Kutu](https://hub.docker.com/r/inductiva/kutu).

We understand that your simulation requirements might extend beyond the 
standard offerings - particularly if you need to run custom code as part 
of a pre- or post-processing workflow.

To give you full flexibility, Inductiva supports **Custom Images**, allowing 
you to run simulations using your own Docker containers.

Simulations are executed using **Apptainer**, an open-source container platform 
built for high-performance computing (HPC). Apptainer uses the 
**Singularity Image Format (SIF)** — a secure, portable, single-file container 
format optimized for HPC workloads.

When running a simulation, Inductiva automatically converts your Docker image into SIF format before executing it. This enables the integration between user-defined Docker environments and our HPC infrastructure.

Your custom images remain **private** — Inductiva securely manages storage and execution of your private containers.

## Create a Singularity Image Format (SIF) File
To use your Docker image with Inductiva, you first need to convert it to the Singularity Image Format (SIF). This can be done easily using the Inductiva CLI:

```
inductiva containers convert <DOCKER IMAGE NAME/ID/URL> <OUTPUT FILE>
```

Or, to reference a remote Docker Hub image directly:

```
inductiva containers convert docker://username/my-private-image:latest private-image.sif
```

### What Happens During Conversion:
- The command locates the specified Docker image in your local Docker registry (or pulls it from Docker Hub if prefixed with `docker://`);
- A conversion container is launched with **Apptainer** pre-installed;
- The Docker image is then converted into a `.sif` file, ready for use in HPC environments.

## Save the SIF format into Inductiva Storage
To use your container in simulations, the `.sif` file must be uploaded to your **Inductiva storage**.

While manual upload is possible, the easiest way is to convert and upload in a single step using the CLI:

```
inductiva containers upload <DOCKER IMAGE NAME/ID/URL>
```

### What This Command Does:
- Convert the Docker image to the `.sif` format;
- Upload the resulting file directly to your private storage on Inductiva.

Once uploaded, your custom container is immediately available for use in simulations via the Inductiva API.

## Run Simulations with a Custom Container
Once your custom Docker image has been uploaded to Inductiva storage in 
`.sif` format, you can reference it directly in your simulation script using the 
`inductiva.simulators.CustomImage` class, as shown below:

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

### Important Notes
- The `inductiva://` prefix is used to reference files stored in your Inductiva private storage;
- The commands list defines what will be executed inside the container;
- You can prepare your input files inside `./input_dir`, which will be mounted inside the container.

Once the simulation is submitted, Inductiva automatically handles the environment setup, container execution, and output management.
