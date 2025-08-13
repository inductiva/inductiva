# Creating a Python Dockerfile for Inductiva

With **Inductiva**, you can run your own software using custom Docker images,
just like in [this example](./integrate-your-docker-container.md).

In this tutorial, we’ll guide you through a few simple steps to run your Python
code on the cloud with Inductiva.

You will start by creating your own **Dockerfile** that builds into a Docker
image containing your Python environment (with all the libraries you need).  
After that, we’ll build the image and send it to Inductiva so you can run your
code in the cloud.

This image will be **available only to you**, and no other users will have access to it.


## Example Dockerfile

```dockerfile
# Start from Ubuntu base image
FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --break-system-packages <all other packages you need>

# Set python3 as the default python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1
```

In this Dockerfile, you will only need to replace `<all other packages you need>` with the
Python packages you want to install.
Example:

```dockerfile
RUN pip3 install --break-system-packages numpy pandas matplotlib
```

## Building the Docker Image

To build the Docker image you just need to run the following command in the same
directory where your Dockerfile is located:

```bash
docker build -t my-python-image --platform linux/amd64 .
```

> **Note**: You can replace `my-python-image` with any name you prefer for your image. The `--platform linux/amd64` flag ensures compatibility with Inductiva's infrastructure, which currently supports only `linux/amd64` images.

## Running the Docker Image on Inductiva
Now that you have your Docker image ready you just need to upload it to Inductiva
with the following command:

```bash
inductiva containers upload my-python-image
```

This command will upload your Docker image to Inductiva, making it available for
you to run simulations using this image.

## Running Your Code
To run your Python code using the Docker image you just uploaded, you can use the
following code snippet:

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c3d-highcpu-16",
    data_disk_gb=20,
    spot=True)

python_image = inductiva.simulators.CustomImage(
    container_image="inductiva://my-containers/my-python-image.sif")


task = python_image.run(
    input_dir="/Users/paulobarbosa/Downloads/empty",
    commands=[
        'python -c "import pandas as pd"',
        'python -c "import numpy as np"'],
    on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()
```

It is this simple! You can now run your Python code on the cloud using your
custom Docker image. The `input_dir` parameter can be set to any directory on your
local machine with your python scripts and data files.

Happy simulations! 🚀

```{banner_small}
:origin: byos-python-docker
```