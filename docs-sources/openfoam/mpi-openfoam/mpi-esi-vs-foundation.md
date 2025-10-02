# MPI Behavior in OpenFOAM-ESI vs OpenFOAM-Foundation

Inductiva supports two OpenFOAM distributions: **ESI** and **Foundation**. While
both achieve the same goal, their internal implementations differ. A notable
difference lies in how they handle the `runParallel` command, specifically in
their use of MPI.

When you run `runParallel`, OpenFOAM uses MPI in the background to execute your
simulation in parallel. However, the way MPI is invoked differs between the two
distributions.

## Behavior by Distribution

### **OpenFOAM-ESI**

* Uses `mpirun` with the `--oversubscribe` flag (we use OpenMPI).
* Relevant script snippet:

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

* Uses `mpirun` without additional flags.
* This means it relies on the default behavior of MPI..
* Relevant script snippet:

```bash
( mpirun -np $nProcs $APP_RUN -parallel "$@" < /dev/null >> log.$LOG_SUFFIX 2>&1 )
```

## What Changes?

The difference in how MPI is invoked leads to distinct behaviors between OpenFOAM-ESI and OpenFOAM-Foundation for the same simulation.

### **OpenFOAM-Foundation**

* Without extra flags, MPI enforces a strict rule: the number of processes requested must be **equal to or less than the number of physical CPU cores** available.
* For example, on a `c2d-highcpu-4` machine (which has **2 physical cores** and **4 vCPUs**), MPI will allow **at most 2 processes**.
* This means your simulation domain can only be split into a number of parts equal to or less than the number of physical cores on your computational resource.

### **OpenFOAM-ESI**

* Uses the `--oversubscribe` flag, which tells MPI to **ignore restrictions on process count**.
* This allows you to run more processes than physical cores, giving you greater flexibility in dividing your simulation domain.
* For example, by default a `c2d-highcpu-4` machine will have 4 vCPUs supported by 2 physical cores. However, OpenFOAM-Foundation will only use at most 2 MPI processes.
* However, running more processes than available vCPUs is **not recommended**, as it can significantly degrade simulation performance.

Keep reading to see how OpenFOAM-Foundation behaves with Inductiva.

```{toctree}
:hidden:
MPI ESI vs Foundation <mpi-esi-vs-foundation>
OpenFOAM-Foundation on Inductiva <foundation-inductiva>
```