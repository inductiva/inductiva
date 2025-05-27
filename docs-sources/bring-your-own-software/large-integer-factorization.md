# Run large integer factorization

All integers can be expressed as a product of prime numbers—a principle known as the fundamental theorem of arithmetic.
Integers with more than one prime factor are called composite numbers, while those with only one prime factor (themselves) are called prime numbers.

Determining the prime factors of a number is known as factorization, a core problem in number theory. Modern cryptographic algorithms, such as RSA, rely on the computational difficulty of factorizing large composite numbers to ensure security.

In this example, we will demonstrate how to factor one of the RSA test numbers using the Inductiva API and the open-source software [CADO-NFS](https://gitlab.inria.fr/cado-nfs/cado-nfs.git), which implements the Number Field Sieve (NFS)—one of the most efficient known algorithms for factoring large integers.

To do this, we will:
	1.	Build a Docker container image with cado-nfs,
	2.	Upload it to Inductiva,
	3.	And execute the factorization on cloud-based compute machines.


## Building the Docker Image

The following Dockerfile installs `cado-nfs` and its dependencies:

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

To build it and save it with the tag `cado-nfs`, run the following command in the directory where the Dockerfile is located:

```bash
docker build -t cado-nfs .
```

If you now run `docker image ls`, you should see the `cado-nfs` image listed.

## Uploading the Docker Image to Inductiva

To upload the Docker image to your Inductiva storage, run the following command:

```bash
inductiva containers upload cado-nfs
```

Note that `cado-nfs` is the name of the image you built in the previous step.
This command will convert the image to Singularity Image Format (SIF) and upload it to your Inductiva storage, making it available for use in simulations.

Run `inductiva containers list` to verify that the image has been uploaded successfully:
```bash
> inductiva containers list                                                                                                                                                          py inductiva lpcunha@lithium

 NAME           SIZE        CREATION TIME
 cado-nfs.sif   337.56 MB   26/05, 10:27:51

```

## Running the Factorization task

Let's check that everything works smoothly by factorizing a number the number `90377629292003121684002147101760858109247336549001090677693` taken from the README of the `cado-nfs` repository.

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

For this first test, we'll use a `cd2-standard-4` machine, which has 4 vCPUs and 16 GB of RAM.
We are using the `CustomImage` "simulator", which allows us to select any Docker image that's publicly available on Docker Hub or that we have uploaded to our Inductiva storage, as we did
above. We are using an `empty` input directory because `cado-nfs` does not require any input files for this task (note that this is actually a directory called empty in my local machine, should we start allowing running tasks without any inputs?).
The `commands` parameter specifies the command to run inside the container, which in this case is the `cado-nfs.py` script followed by the number we want to factorize.
To keep things tidy, we are also specifying a project name `cado-nfs` to group this and related tasks for easier management and retrieval later on.

`task.wait()` will block until the task is completed, showing the command output in real-time, and `task.print_summary()` will print a summary of the task:

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

Running the factorization took approximately 21 seconds, excluding time spent setting up the machine and preparing the environment.

In the `stdout.txt` file, that you can find in the console or by downloading the task outputs via the CLI (`inductiva tasks download <task_id>`), you will find the factorization result:

```
260938498861057 588120598053661 760926063870977 773951836515617
````

So, `90377629292003121684002147101760858109247336549001090677693` can be factorized into the product of these four prime numbers.


## Factorizing examples from the RSA Challenge

Now that we've seen that everything is running smoothly, let's bump up the challenge, and try to factor some of the numbers from the [RSA Challenge](https://en.wikipedia.org/wiki/RSA_numbers).
These numbers are semiprime, meaning they can be expressed as the product of two, also very large, prime numbers.

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

We compiled a list of RSA numbers to factor, which you can find here. Specifically, we selected RSA-100, RSA-110, …, RSA-150, where the number after the dash indicates the number of decimal digits in the RSA number.

To perform the factorization, we created an `ElasticMachineGroup` using the c2d-highcpu-112 machine type, which provides 112 vCPUs and 224 GB of RAM—sufficient resources to handle the computational demands of these tasks.

The machine group is configured with a maximum number of machines equal to the number of RSA numbers, ensuring that each number can be processed in parallel on a separate machine.

Thanks to its elastic nature, the group automatically scales up and down: a new machine is provisioned for each task when it starts, and terminated as soon as the task finishes. This means machines do not wait for all tasks to complete—they run independently and shut down immediately after completing their assigned job.

Here are some statistics of the tasks we ran:


|   RSA Number  |     Execution Time     |   Estimated Cost   |
|:-------------:|:----------------------:|:--------:|
|  RSA-100  |  4 minutes and 43 seconds | 0.057 US$ |
|  RSA-110 | 8 minutes and 6 seconds | 0.098 US$ |
|  RSA-120 | 16 minutes and 39 seconds   | 0.19 US$ |
|  RSA-130 |   42 minutes and 1 second  | 0.50 US$ |
|  RSA-140 |  2 hours and 12 minutes    | 1.58 US$ |
|  RSA-150 |    ----  | --- |



## Conclusion

In this tutorial, we demonstrated how one can leverage the Inductiva API to easily run arbitrary computational tasks in the Cloud,
using as an example the factorization of integers using the CADO-NFS software.









