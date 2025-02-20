---
orphan: true
---

# Preparation

Before running a simulation in OpenFAST, it’s crucial to ensure that all
necessary files are properly set up. Here we will guide you through the preparation
process. We will walk you through downloading the required input files, and
generating any additional components needed, such as the `DISCON_OC3Hywind.dll`
file. By the end of this section, you’ll have everything in place to proceed
with running the simulation.

Let’s get started!  


## Our use case: 5MW_OC4Semi_WSt_WavesWN

In this tutorial, we will show you how to do this using the 
"5MW_OC4Semi_WSt_WavesWN" example, available from the OpenFast 
[GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.
 
This example is an extension of the reference case described in 
["Definition of a 5-MW Reference Wind Turbine for Offshore
System Development"](https://www.nrel.gov/docs/fy09osti/38060.pdf).

All files needed are available from OpenFast GitHub so
let's get started.

### Requirements: Setting up your files

Before we start running the simulation, we need to ensure that all required
files are properly set up. In this section, we'll go through the steps to
download, and prepare the necessary input files.

#### Step 1: Downloading you simulation files

We are going to run the `5MW_OC4Semi_WSt_WavesWN` case as
it is originally defined in the [GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.

We start by downloading the folders `5MW_OC4Semi_WSt_WavesWN` and the `5MW_Baseline` to our local input directory. This folders are located [here](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast).

After downloading the folders and moving them your `input_files` should look something like this:

- input_files
    - 5MW_Baseline
    - 5MW_OC4Semi_WSt_WavesWN

#### Step 2: Building `DISCON_OC3Hywind.dll`

After downloading the files we need to build the `DISCON_OC3Hywind.dll` file and place it in `5MW_Baseline/ServoData/`. 

To do so we need to do the following steps:
```
cd input_files/5MW_Baseline/ServoData/DISCON_OC3/
mkdir build
cd build
cmake ..
make
mv DISCON_OC3Hywind.dll ../..
```

We should now have our input directory looking like this:
  
- input_files
    - 5MW_Baseline
        - ServoData
            - DISCON_OC3Hywind.dll
    - 5MW_OC4Semi_WSt_WavesWN

> Note: you can also download de dll file [here](https://storage.googleapis.com/inductiva-simulators-sources/DISCON_OC3Hywind.dll). 
Don't forget to paste this file in the 5MW_Baseline/ServoData folder.

You now have all the necessary files to run your simulation.