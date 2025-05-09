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
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
