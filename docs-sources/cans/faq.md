**Find answers to commonly asked questions about CaNS.**

<br>

# FAQ

## 1. What should I do if I encounter the "NCCL error (error code 5)" and "Invalid usage (invalid grid descriptor)" when running CaNS with GPU support?
If you encounter the following errors when running CaNS with GPU support:

```
CUDECOMP:ERROR: /CaNS/dependencies/cuDecomp/src/cudecomp.cc:69 NCCL error. (error code 5)
CUDECOMP:ERROR: /CaNS/dependencies/cuDecomp/src/cudecomp.cc:147 Invalid usage. (invalid grid descriptor)
```

You can resolve the issue by changing the settings in your input files.

In the configuration file, locate the following lines:

```
cudecomp_is_t_enable_nccl = T
cudecomp_is_h_enable_nccl = T
```

Then, change them to:

```
cudecomp_is_t_enable_nccl = F
cudecomp_is_h_enable_nccl = F
```

This should eliminate the error and allow the simulation to run smoothly.

<br>

## 2. I want to run my simulation on a GPU-capable machine. How many `n_vcpus` should I use?

CaNS can be sensitive when running on GPUs, so it's important to set the `n_vcpus`
parameter correctly to avoid issues. As a general rule, set `n_vcpus` equal to
the total number of GPUs being used for the simulation.

**Examples:**

* If you're running on a single machine with **2 GPUs**, set `n_vcpus = 2`.
* If you're running on a cluster with **4 machines**, each with **1 GPU**, set `n_vcpus = 4`.

This ensures that each GPU is correctly assigned and utilized during the simulation.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
