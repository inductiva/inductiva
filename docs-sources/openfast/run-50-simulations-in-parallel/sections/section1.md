# Prerequisites
Before running an OpenFAST simulation, it's important to ensure that all the necessary files are set up correctly. 
Here we will guide you through the preparation process. We will guide you through downloading the required input files and generating any additional components, 
such as the `DISCON_OC3Hywind.dll` file. By the end of this section you'll have everything you need to run the simulation.

Let’s get started!

### Step 1: Download your simulation files
Download the required [`5MW_OC4Semi_WSt_WavesWN`](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) and [`5MW_Baseline`](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_Baseline) folders and place them in your working directory. Your `input_files` should now look like this:
- input_files
    - 5MW_Baseline
    - 5MW_OC4Semi_WSt_WavesWN

### Step 2: Build the DLL file
Now we need to build the `DISCON_OC3Hywind.dll` file and place it in `5MW_Baseline/ServoData/`.

**Want to skip this step?** Download the DLL file [here](https://storage.googleapis.com/inductiva-simulators-sources/DISCON_OC3Hywind.dll). 
Don't forget to place it in the `5MW_Baseline/ServoData` folder.

If you wish to build the file yourself, please follow the steps outlined in this [tutorial](https://inductiva.ai/guides/openfast/build-dll-file).

Your input directory should look like this:
- input_files
    - 5MW_Baseline
        - ServoData
            - DISCON_OC3Hywind.dll
    - 5MW_OC4Semi_WSt_WavesWN

You're ready to send your simulation to the Cloud!