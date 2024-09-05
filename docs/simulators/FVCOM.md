# FVCOM

[FVCOM](https://github.com/CaNS-World/CaNS) (Finite Volume Community Ocean Model)
is a three-dimensional hydrodynamic model designed for simulating coastal and
ocean dynamics. It uses an unstructured grid with finite-volume methods, allowing
for greater flexibility in representing complex coastlines, estuaries, and variable
bathymetry. FVCOM is particularly effective in modeling ocean circulation, tides,
and coastal processes, providing high-resolution simulations of water currents,
temperature, salinity, and ecosystem interactions. Its versatility and accuracy
make it a valuable tool for studying coastal environments, marine habitats, and
climate change impacts on ocean systems.

FVCOM has been compiled using OpenMPI 4.1.2 and with the following flags:
 - FLAG_USE_NETCDF4     (DUSE_NETCDF)
 - FLAG_1               (DDOUBLE_PRECISION)
 - FLAG_4               (DMULTIPROCESSOR)
 - FLAG_411             (DMETIS_5)
 - FLAG_8               (DLIMITED_NO)
 - FLAG_10              (DGCN)
 - FLAG_44              (DTVD)
 - FLAG_15              (DMPDATA)

Any flags not mentioned above are currently disabled. This is our only compiled
version of FVCOM at the moment, but we are actively working on a solution that
will allow users to customize and specify their own flags for on-the-fly compilation.

## Example

```python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fvcom-input-example.zip", unzip=True)

# Initialize the Simulator
fvcom = inductiva.simulators.FVCOM()

# Run simulation with config files in the input directory
task = fvcom.run( input_dir=input_dir,
                  working_dir="run/",
                  case_name="tst",
                  debug=7,
                  n_vcpus=1)

task.wait()
task.download_outputs()

```

**Closing Notes**: There is currently a bug affecting the timezone argument in
the `nml` file. If you encounter issues, please set the timezone to either `None`
or `UTC` in the `namelist.nml` file as a workaround.
[More information](https://github.com/FVCOM-GitHub/FVCOM/issues/27)

## What to read next

If you are interested in FVCOM, you may also be interested in checking
the following related simulators that are also avaiable via Inductiva API:

* [SWAN](SWAN.md)
* [SWASH](SWASH.md)
* [OpenFAST](OpenFAST.md)
