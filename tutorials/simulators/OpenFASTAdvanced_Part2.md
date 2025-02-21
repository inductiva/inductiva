---
orphan: true
---

# Preparation

Before running a simulation in OpenFAST, it's crucial to ensure that all
necessary files are properly set up. Here we will guide you through the preparation
process. We will walk you through downloading the required input files, and
generating any additional components needed, such as the `DISCON_OC3Hywind.dll`
file. By the end of this section, you'll have everything in place to proceed
with running the simulation.

Let's get started!

## Preparation Step 1: Downloading your simulation files

We are going to run the `5MW_OC4Semi_WSt_WavesWN` case as
it is originally defined in the [GitHub](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) page.

We start by downloading the folders `5MW_OC4Semi_WSt_WavesWN` and the `5MW_Baseline` to our local input directory. This folders are located [here](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast).

After downloading the folders and moving them your `input_files` should look something like this:

- input_files
    - 5MW_Baseline
    - 5MW_OC4Semi_WSt_WavesWN

## Preparation Step 2: Building `DISCON_OC3Hywind.dll`

After downloading the files we need to build the `DISCON_OC3Hywind.dll` file and
place it in `5MW_Baseline/ServoData/`.

> Note: You can also download de dll file [here](https://storage.googleapis.com/inductiva-simulators-sources/DISCON_OC3Hywind.dll). By downloading the file you can skip the following sections.
Don't forget to paste this file in the 5MW_Baseline/ServoData folder.

### Prerequisites

Before proceeding, ensure you have the following installed:

- A Fortran compiler (`gfortran` can be installed with `sudo apt-get install gfortran`)
- `make` (included in `build-essential`, install with `sudo apt-get install build-essential`)

> Note: This commands are intended for Ubuntu. If you're using a different
operating system, you may need to adjust them accordingly.

### Building the DLL

Once the prerequisites are met, run the following commands to generate the
`DISCON_OC3Hywind.dll` file:
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

You now have all the necessary files to run your simulation.

[Running the simulation](OpenFASTAdvanced_Part3.md)