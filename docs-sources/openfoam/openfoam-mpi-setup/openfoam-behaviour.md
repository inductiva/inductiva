# MPI Behaviour in OpenFOAM
Inductiva supports two OpenFOAM distributions: **ESI** and **Foundation**. Both use MPI in the background when executing the `runParallel` command for parallel simulations, but they differ in how MPI is invoked.

## Behaviour by Distribution

### OpenFOAM-ESI
Uses `mpirun` with the `--oversubscribe` flag (OpenMPI is used). Take a look at how this is handled in the `runParallel` script:

```bash
local mpirun="mpirun"
case "$FOAM_MPI" in
(msmpi*)
    mpirun="mpiexec"
    ;;
(*openmpi*)
    mpiopts="--oversubscribe"
    ;;
esac
```

### **OpenFOAM-Foundation**
Uses `mpirun` without additional flags, relying on MPI’s default behavior — launching one process per physical core (not per available thread) Take a look at how this is handled in the `runParallel` script:

```bash
( mpirun -np $nProcs $APP_RUN -parallel "$@" < /dev/null >> log.$LOG_SUFFIX 2>&1 )
```

## Impact of the Differences
The variation in how MPI is invoked causes different behaviours between OpenFOAM-ESI and OpenFOAM-Foundation when running the same simulation. The following sections explore these differences in more detail.

> Learn more about MPI on computational resources [here](https://inductiva.ai/guides/how-it-works/machines/mpi-on-vms).

```{banner_small}
:origin: openfoam
```

```{toctree}
:hidden:
sections/openfoam-esi.md
sections/openfoam-foundation.md
```
