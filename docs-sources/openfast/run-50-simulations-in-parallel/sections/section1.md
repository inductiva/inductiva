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

### Step 2: Specify the Correct DLL File

This simulation requires the `DISCON_OC3Hywind.dll` file. Normally, this DLL
must be compiled from the source code provided in the `5MW_Baseline/ServoData` folder.

To simplify this process, we’ve pre-compiled all three commonly used DLLs and bundled them directly in the OpenFAST Docker images. You can find them at:

* `/DLLs/DISCON.dll`
* `/DLLs/DISCON_ITIBarge.dll`
* `/DLLs/DISCON_OC3Hywind.dll`

Once you know the location of the correct DLL, simply update your `NRELOffshrBsline5MW_OC4DeepCwindSemi_ServoDyn.dat` file to point to it:

```diff
---------------------- BLADED INTERFACE ---------------------------------------- [used only with Bladed Interface]
-"../5MW_Baseline/ServoData/DISCON_OC3Hywind.dll"    DLL_FileName
+"/DLLs/DISCON_OC3Hywind.dll"    DLL_FileName
```

With this small adjustment, your simulation is ready to run in the cloud!
