# Prerequisites

## 1. Download the Simulation Data
Download the simulation from [NASA's collection bundle](https://data.ornldaac.earthdata.nasa.gov/protected/bundle/DeltaX_Delft3D_294_Terrebonne_2303.zip) *(registration required)*.

After unzipping, the folder structure should look like this:

```
DeltaX_Delft3D_294_Terrebonne_2303
   ├─ comp
   ├─ data
   └─ guide
```

## 2. Select the Required Files
The configuration files are located in the `data/` directory. Several NetCDF output files (`.nc4`) are included but **not needed** for this tutorial:

```
data
  ├─ Fall2021_Delft3D_setup_294.zip
  ├─ Spring2021_Delft3D_setup_294.zip
  ├─ site294_IMAR.nc4
  ├─ site294_d3d_output_ssc_mud_Fall2021.nc4
  └─ site294_d3d_output_wl_Spring2021.nc4
```

Each zip archive contains a complete Delft3D model setup, one for the Spring and one for the Fall deployment.

Unzip both archives and organize them into a new folder named `my_project`, placed inside the main `DeltaX_Delft3D_294_Terrebonne_2303/` directory:

```
DeltaX_Delft3D_294_Terrebonne_2303
   ├─ comp
   ├─ data
   ├─ guide
   └─ my_project
        ├─ Fall2021_Delft3D_setup
        └─ Spring2021_Delft3D_setup

```

## 4. Add the `config_d_hydro.xml`
Delft3D requires an XML configuration file to trigger the simulation. This file is **not included** in the dataset, so you will need to add one. You can use the following template:

```
<?xml version="1.0" encoding="iso-8859-1"?>
<deltaresHydro xmlns="http://schemas.deltares.nl/deltaresHydro" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://schemas.deltares.nl/deltaresHydro http://content.oss.deltares.nl/schemas/d_hydro-1.00.xsd">
    <flow2D3D name="inductiva_flow">
        <library>flow2d3d</library>
        <mdfFile>site_294.mdf</mdfFile>
    </flow2D3D>
    <delftOnline>
        <enabled>false</enabled>
    </delftOnline>
</deltaresHydro>
```

Save this file as `config_d_hydro.xml` and copy it into both the `Fall2021_Delft3D_setup` and `Spring2021_Delft3D_setup` folders inside your `my_project/`. This ensures that both simulations can be executed.

## Folder Overview
After adding the configuration files, your directory structure should look like this:

```
DeltaX_Delft3D_294_Terrebonne_2303
   ├─ comp
   ├─ data
   ├─ guide
   └─ my_project
        ├─ Fall2021_Delft3D_setup
        │    ├─ Delft3D input files...
        │    └─ config_d_hydro.xml
        │
        └─ Spring2021_Delft3D_setup
             ├─ Delft3D input files...
             └─ config_d_hydro.xml
```

All required data and configuration files are now in place. Next, we’ll prepare the Python script to run the 2021 Spring deployment.