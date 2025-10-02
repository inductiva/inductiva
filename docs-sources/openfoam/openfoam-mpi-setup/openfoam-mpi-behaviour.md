# MPI Behaviour in OpenFOAM
Inductiva supports two OpenFOAM distributions: **ESI** and **Foundation**. While
both achieve the same goal, their internal implementations differ. The difference 
lies in how they handle the `runParallel` command, specifically in their use of MPI.

When executing `runParallel`, OpenFOAM uses MPI in the background to run simulation in parallel. 
However, the method of invoking MPI varies between these two distributions.

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
Uses `mpirun` without additional flags, which means it relies on MPI’s default behaviour. Take a look at how this is handled in the `runParallel` script:

```bash
( mpirun -np $nProcs $APP_RUN -parallel "$@" < /dev/null >> log.$LOG_SUFFIX 2>&1 )
```

## Impact of the Differences
The variation in how MPI is invoked causes different behaviours between OpenFOAM-ESI and OpenFOAM-Foundation when running the same simulation. The following sections explore these differences in more detail.

```{banner_small}
:origin: openfoam
```

```{toctree}
:hidden:
sections/openfoam-esi.md
sections/openfoam-foundation.md
```
