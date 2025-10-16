# Prerequisites

## 1. Download the Simulation Data
Download the simulation data from [Zenodo](https://zenodo.org/records/8308164). The files are compatible with **SWAN version 41.45**.

After unzipping, the folder structure should look like this:

```
SWAN_TerraNovaBayPolynya_files/
  ├── Tp_from_satellite_swangrid/
  └── swan/
       ├── input/
       ├── S0/
       └── S2_f5/
```

## 2. Select the Required Files
For this tutorial, you only need the following:
- The `input/` folder (wind, ice concentration, and grid data)
- `.swn` input files from `S0/` and `S2_f5/` (one per polynya event)

The `S0/` and `S2_f5/` folders also contain large output files (`.mat`, `.dat`, `.sp1`, `.sp2`). **Do not copy these**. Only the `.swn` files are required for each simulation.

The `Tp_from_satellite_swangrid` folder contains satellite-derived peak periods and directions, interpolated onto the SWAN grid. These were used in the original study to validate and calibrate the model (bias, RMSD, spectral tail comparisons), but are **not required** for this tutorial.

## 3. Organize the `tutorial/` Folder
Create a new directory named `tutorial/` and copy only the required files into it. Your structure should look like this:

```
tutorial
   |
   + polynya
      |
      + S0
      |  |
      |  + polynya2D_20161005.swn
      |  |	
      |  + polynya2D_20161024.swn
      |  |
      |  ...
      |  |
      |  + polynya2D_20201019.swn
      |
      |
      + S2_f5
      |  |
      |  + polynya2D_20161005.swn
      |  |	
      |  + polynya2D_20161024.swn
      |  |
      |  ...
      |  |
      |  + polynya2D_20201019.swn
      |
      + input
         |
         + Aice0_20161005_dx200.dat 
         |
         ...
         |
         + wnd_20211007_10m.dat
```

## 4. Modify Each `.swn` File
Before running any scripts, you need to edit each `.swn` input file in both the `S0/` and `S2_f5/` folders to prevent a known Fortran runtime issue.

In every `.swn` file, locate the following line:

```
SPECOUT 'COMPGRID' SPEC1D ABS '20211007.sp1'
```

and **comment it out** by adding an exclamation mark (!) at the beginning:

```
!SPECOUT 'COMPGRID' SPEC1D ABS '20211007.sp1'
```

This modification avoids a common Fortran runtime error that can occur when SWAN attempts to output large spectral files over the full computational grid.

## Folder Overview
Here's what each subfolder in `tutorial/polynya` contains:
- `swan/input`: Grid definitions, wind fields, and ice concentration maps for 10 polynya events. These are the raw fields for your SWAN simulations.
- `swan/S0`: Input files for the open water baseline scenario (S0). These runs include only wind-forced waves, with no ice effects.
- `swan/S2_f5`: Input files for the full ice-influenced scenario (S2_f5), which included frazil/grease effects and applies ice and a power-law exponent of 5 are applied.

All required data is now in place. Next, let's prepare the Python script to run the baseline S0 scenario.

---
**Notas para o Paulo**:

**Verificar que erro é este, para perceber se isto tem de ser feito**:

"
Before actually starting the script, we need to change one line in the 
S0/polynya2D_20211007.swn to avoid a problem explained (here). The following line:

```
SPECOUT 'COMPGRID' SPEC1D ABS '20211007.sp1'
```

needs to be commented out:

```
!SPECOUT 'COMPGRID' SPEC1D ABS '20211007.sp1'
```
to avoid a common Fortran runtime error that may happen when SWAN tries to
output large spectrum files over the entire computational grid. 
"

**O erro encontra-se abaixo, incluído em FAQ. Valida isto.**

"
Q: Why am I getting a "Fortran runtime error: End of file", at the very end of a
simulator that seemed to be going ok?

A: If you are getting a Failed status at the end of what seems to be an otherwise 
successful simulation, you may be encountering a common problem that happens when 
running SWAN in parallel / multi-threaded setting that produces large 1D/2D 
spectrum files.

Check in the logs to see if you find an error message such as:

```
At line 4880 of file /tmp/swan/src/swanparll.f (unit = 698, file = 'out_spectrum.sp1-086')
Fortran runtime error: End of file

Error termination. Backtrace:
#0  0x7f3a8b575960 in ???
#1  0x7f3a8b5764d9 in ???
#2  0x7f3a8b7ca17b in ???
#3  0x7f3a8b7ca752 in ???
#4  0x7f3a8b7c6f1b in ???
#5  0x7f3a8b7cbe3c in ???
#6  0x7f3a8b7cce55 in ???
#7  0x557be791ed88 in swcolspc_
#8  0x557be792560f in swcolout_
#9  0x557be7847aed in swmain_
#10  0x557be7847ebf in main
--------------------------------------------------------------------------
Primary job  terminated normally, but 1 process returned
a non-zero exit code. Per user-direction, the job has been aborted.
--------------------------------------------------------------------------
```

This is a common problem when using SPEC1D/SPEC2D on COMPGRID for an entire 
relatively large grid (e.g. 400×300) that will produce large data and many
writes. If SWAN is running over multiple MPI processes then the volume of writes 
over the many output streams (one per MPI process) can increases the chance of 
I/O contention or hitting file-system limits.

You can overcome this issue by reducing the level of parallelism (e.g. running 
SWAN over smaller number of MPI processes) or, if you need spectra only at a 
handful of points,by defining a small set of POINTS or LINES for output rather 
than the whole COMPGRID.
"
