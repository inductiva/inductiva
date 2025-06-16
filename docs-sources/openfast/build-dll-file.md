# Build a DLL file
In the [Run 50 Simulations in Parallel](https://inductiva.ai/guides/openfast/OpenFAST_advanced) tutorial, we used the `5MW_OC4Semi_WSt_WavesWN` example, an extension of the reference case from the “Definition of a 5-MW Reference Wind Turbine for Offshore System Development”, which can be found on the [OpenFAST GitHub repository](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN).

In some cases, including the **OC3 Hywind**, it is necessary to build the `DISCON_OC3Hywind.dll` file.

To do so, ensure you have the following installed:
- A Fortran compiler (`gfortran` can be installed with `sudo apt-get install gfortran`)
- The `make` software (included in `build-essential`, install with `sudo apt-get install build-essential`)

> Note: These commands are for Ubuntu. If you're using a different operating system, you may need to adapt them.

Now, run the following commands to generate the `DISCON_OC3Hywind.dll` file:
```
cd input_files/5MW_Baseline/ServoData/DISCON_OC3/
mkdir build
cd build
cmake ..
make
mv DISCON_OC3Hywind.dll ../..
```

That's it! You're all set.


