# OpenFOAM-ESI
OpenFOAM-ESI uses the `--oversubscribe` flag with `mpirun`, which instructs MPI to **ignore restrictions on the number of processes** relative to available CPU cores.

This allows you to run more processes than physical cores, which can be helpful for enabling finer domain decompositions.

For example, a `c2d-highcpu-16` machine has 16 vCPUs supported by 8 physical cores. Thanks to `--oversubscribe`, OpenFOAM-ESI lets you run up to 8 MPI processes on this machine.

> ⚠️ Note: While this flexibility can be useful, running more processes than available vCPUs is generally not recommended, as it may significantly degrade simulation performance.