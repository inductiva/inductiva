# Integrate Your Docker Container
One of Inductiva's standout features is its ability to seamlessly run any Docker container directly in the cloud. In this guide, we’ll walk you through the integration process using [FFMPEG](https://ffmpeg.org) as an example. This approach can also be applied to your **own custom Docker images**.

## What You'll Learn
By the end of this guide, you’ll know how to:
- Use a standard Docker image (FFMPEG)
- Upload it to Inductiva’s container registry
- Run it in the cloud with scalable compute

Let's dive in!

## Step 1: Pull the Docker Image
You can use any Docker image from Docker Hub or your own private registry. For this guide, we will use a pre-built FFMPEG image, a widely-used tool for media processing. To get started, run the following command:

```
docker pull linuxserver/ffmpeg --platform linux/amd64
```

## Step 2: Upload to Inductiva
Uploading Docker images to Inductiva is simple — just **one command**. When you upload your image, it will automatically be converted into the [Singularity Image Format](https://en.wikipedia.org/wiki/Singularity_(software)) (.sif), which is commonly used in high-performance computing for its portability and security.

To upload your image, run:

```bash
inductiva containers upload linuxserver/ffmpeg
```

> 💡 Tip: Use `inductiva containers list` to verify the containers you’ve uploaded.

## Step 3: Run It in the Cloud
Once your container is uploaded, you can use it for cloud-based tasks just like any native Inductiva simulator.

Let’s say you want to run a video processing command that reverses and negates the frames of a video while also reversing the audio. Here’s the command:

```bash
ffmpeg -i input.mp4 -vf 'reverse, negate' -af areverse output.mp4
```

Copy and execute the following Python script to run this remotely on the cloud using the Inductiva API:

```python
import inductiva

# Define the machine group
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2d-highcpu-8",
    spot=True,
)

# Define your simulator using the uploaded container
ffmpeg = inductiva.simulators.CustomImage(
    container_image="inductiva://my-containers/ffmpeg.sif"
)

# Run the task
task = ffmpeg.run(
    input_dir="input",  # folder containing 'input.mp4'
    commands=[
        "ffmpeg -i input.mp4 -vf 'reverse, negate' -af areverse output.mp4"
    ],
    on=machine_group,
)

# Wait for it to complete and download the result
task.wait()
task.download_outputs()

# Clean up your resources
machine_group.terminate()
```

> **Note**: Ensure that the `input_dir` includes your `input.mp4` file and the processed video will be saved as output.mp4.

## Why Choose Inductiva for Your Docker Workflows
Inductiva streamlines the process of running Docker containers in the cloud, and here’s why it’s the ideal solution::
- **No changes required to your application**: Just upload your Docker image and you're ready to go
- **No manual VM setup**: Inductiva handles everything automatically.
- **Run any container**: From media processing tools like FFMPEG to PyTorch or custom AI tools.

Whether you’re processing videos, running simulations, or deploying research tools, Inductiva simplifies the process of scaling your Docker workflows to the cloud — without the hassle.

