# Run Large Integer Factorization
Every integer can be expressed as a product of prime numbers — a principle known as the Fundamental
Theorem of Arithmetic.

Integers with more than one prime factor are called composite numbers, while those with only themselves as
a prime factor are called prime numbers.

Determining the prime factors of a number is known as **factorization**, a core problem in number theory.
Modern cryptographic algorithms, such as RSA, rely on the computational difficulty of factorizing large
composite numbers to ensure security.

In this guide, we will demonstrate how to factor one of the RSA test numbers using the Inductiva API
and the open-source software [CADO-NFS](https://gitlab.inria.fr/cado-nfs/cado-nfs.git), which implements
the Number Field Sieve (NFS) — one of the most efficient known algorithms for factoring large integers.

We’ll walk through:
1.	Building a Docker container image with CADO-NFS
2.	Uploading it to Inductiva
3.	Executing the factorization on cloud machines

## Building the Docker Image
First, create a Dockerfile that installs `cado-nfs` and its dependencies:

```dockerfile
FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    cmake \
    libgmp-dev \
    python3 \
    python3-flask \
    python3-requests

RUN git clone https://gitlab.inria.fr/cado-nfs/cado-nfs.git /cado-nfs
WORKDIR /cado-nfs

ENV PREFIX=/usr/local

RUN make cmake
RUN make -j
RUN make install
```

To build the image and save it with the tag `cado-nfs`, run the following command in the directory
where the Dockerfile is located:

```bash
docker build -t cado-nfs .
```

You can verify the image was created by running `docker image ls`. The `cado-nfs` image will be listed.

## Uploading the Docker Image to Inductiva
Upload the Docker image to Inductiva storage by running the following command:

```bash
inductiva containers upload cado-nfs
```

This command will convert the image to Singularity Image Format (SIF) and upload it to your Inductiva storage,
making it available for use in simulations.

Note that `cado-nfs` is the name of the image you built in the previous step.

Run `inductiva containers list` to verify that the image has been uploaded successfully:

```bash
> inductiva containers list

 NAME           SIZE        CREATION TIME
 cado-nfs.sif   337.56 MB   26/05, 10:27:51

```

## Running the Factorization Task
We’ll now factor a sample number from the `cado-nfs` repository:

```python
import inductiva

cloud_machine = inductiva.resources.MachineGroup("c2d-highcpu-4")

cado_nfs = inductiva.simulators.CustomImage(
    "inductiva://my-containers/cado-nfs.sif")

task = cado_nfs.run(
    on=cloud_machine,
    input_dir="empty",
    commands=[
        "cado-nfs.py 90377629292003121684002147101760858109247336549001090677693"
    ],
    project="cado-nfs",
)
task.wait()
task.print_summary()
```

For this initial test, we use a `cd2-standard-4` machine, which has 4 vCPUs and 16 GB of RAM.

We are using the `CustomImage` simulator, which lets us select any Docker image, either one publicly available on Docker Hub or one we’ve uploaded to our Inductiva storage, as demonstrated earlier.

In this case, we are using an `empty` input directory because `cado-nfs` does not require any input files for this task.

The `commands` parameter specifies the command to run inside the container. Here, it's the `cado-nfs.py` script followed by the number we want to factor.

To keep things organized, we’ve assigned the project name `cado-nfs` to this task, making it easier to manage and retrieve related runs later on.

Calling `task.wait()` will block execution until the task completes, showing the command's output in real time. Once done, task.print_summary() will display a summary of the task, as follows:

```
Task status: Success

Timeline:
	Waiting for Input         at 26/05, 17:43:55      0.692 s
	In Queue                  at 26/05, 17:43:56      35.529 s
	Preparing to Compute      at 26/05, 17:44:32      2.063 s
	In Progress               at 26/05, 17:44:34      21.255 s
		└> 21.138 s        cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
	Finalizing                at 26/05, 17:44:55      0.348 s
	Success                   at 26/05, 17:44:55

Data:
	Size of zipped output:    5.06 KB
	Size of unzipped output:  37.56 KB
	Number of output files:   2

Estimated computation cost (US$): 0.00024 US$
```

As you can see in the “In Progress” line, the part of the timeline that represents the actual execution of the simulation, the core computation time of this simulation was approximately 21 seconds.

The factorization result can be found in the `stdout.txt` file, which is accessible through the Inductiva Console or by downloading the task outputs using the CLI (`inductiva tasks download <task_id>`):

```
260938498861057 588120598053661 760926063870977 773951836515617
```

So, `90377629292003121684002147101760858109247336549001090677693` can be factorized into the product of these four prime numbers.


## Factorizing Examples from the RSA Challenge
Now let's bump up the challenge and attempt to factor some of the numbers from the [RSA Challenge](https://en.wikipedia.org/wiki/RSA_numbers). These numbers are semiprime, meaning each is the product of two (very large) prime numbers.

```python
import inductiva

rsa_numbers = [
    "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139",
    "35794234179725868774991807832568455403003778024228226193532908190484670252364677411513516111204504060317568667",
    "227010481295437363334259960947493668895875336466084780038173258247009162675779735389791151574049166747880487470296548479",
    "1807082088687404805951656164405905566278102516769401349170127021450056662540244048387341127590812303371781887966563182013214880557",
    "21290246318258757547497882016271517497806703963277216278233383215381949984056495911366573853021918316783107387995317230889569230873441936471",
    "155089812478348440509606754370011861770654545830995430655466945774312632703463465954363335027577729025391453996787414027003501631772186840890795964683",
]

cloud_machine = inductiva.resources.ElasticMachineGroup(
    "c2d-highcpu-112",
    max_machines=len(rsa_numbers),
)
cado_nfs = inductiva.simulators.CustomImage(
    "inductiva://my-containers/cado-nfs.sif")

for rsa_number in rsa_numbers:
    task = cado_nfs.run(
        on=cloud_machine,
        input_dir="empty",
        commands=[f"cado-nfs.py {rsa_number}"],
        project="cado-nfs-2",
    )
    task.set_metadata({
        "rsa_number": str(len(rsa_number)),
    })
```

We compiled a list of RSA numbers to factor, specifically RSA-100, RSA-110, ..., RSA-150. The number after the dash indicates the number of decimal digits in the corresponding RSA number.

To tackle the computation, we created an `ElasticMachineGroup` using the `c2d-highcpu-112` machine type, which provides 112 vCPUs and 224 GB of RAM—sufficient resources for these intensive tasks.

The machine group is configured to allow up to as many machines as there are RSA numbers, ensuring that each number can be processed in parallel on a separate machine.

Thanks to the elastic nature of this setup, a new machine is provisioned as soon as a task starts and is automatically shut down upon completion. This ensures efficient resource usage: machines don’t idle while waiting for others to complete their assigned job.

Here’s a summary of the performance and cost of the RSA factorization tasks:

|   RSA Number  |     Execution Time     |   Estimated Cost   |
|:-------------:|:----------------------:|:--------:|
|  RSA-100  |  4 minutes and 43 seconds | 0.057 US$ |
|  RSA-110 | 8 minutes and 6 seconds | 0.098 US$ |
|  RSA-120 | 16 minutes and 39 seconds   | 0.19 US$ |
|  RSA-130 |   42 minutes and 1 second  | 0.50 US$ |
|  RSA-140 |  2 hours and 12 minutes    | 1.58 US$ |
|  RSA-150 |    6 hours and 39 minutes  | 4.87 US$ |

## Wrapping Up
In this guide, you learned how to use the Inductiva API to run arbitrary computational tasks in the cloud. As a concrete example, we used the CADO-NFS software to factor large semiprime integers, including several from the RSA Challenge.

This demonstrates how easily custom workloads can be executed with Inductiva, unlocking access to powerful computational resources with minimal setup.







