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

In each `.swn` file, look for a line that starts with:

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