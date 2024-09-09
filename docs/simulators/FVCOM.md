# FVCOM

[FVCOM](https://www.fvcom.org/) (Finite Volume Community Ocean Model)
is a 3D hydrodynamic model designed for simulating coastal and ocean dynamics.
It utilizes an unstructured grid with finite-volume methods, providing flexibility
in modeling complex coastlines, estuaries, and varying bathymetry. FVCOM excels
in simulating ocean circulation, tides, and coastal processes, offering
high-resolution outputs for water currents, temperature, salinity, and ecosystem
interactions. Its versatility and precision make it an essential tool for
studying coastal environments, marine ecosystems, and the impacts of climate
change on ocean systems.

We have compiled two versions of FVCOM: the standard **fvcom** binary and an
additional **fvcom_estuary** binary, which is configured to run the Estuary test
case included in the FVCOM package. Users can check the compilation modules used
in each binary in the respective `make.inc` files located at `/make.inc` and
`/make_estuary.inc`.

These are the only FVCOM versions available at the moment, but we are working on
a solution to allow users to customize and specify their own flags for on-the-fly
compilation.

If you are having trouble with the input namelist file, you can do the following to
generate a valid `nml` file in the `working_dir` with the expected format for
the requested model:

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
                  create_namelist="tst",
                  n_vcpus=1)

task.wait()
task.download_outputs()

```

## Example

In the following example, we run the default model with a debug level of 7 and 1
MPI process on a very small test scenario. This is just a simple example to check
if the simulator is working correctly.

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
