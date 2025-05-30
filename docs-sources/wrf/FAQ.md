**Find answers to commonly asked questions about WRF.**

<br>

# FAQ

## In your tutorials, you're using the case `em_real`. What other cases are available?

WRF includes several pre-configured test cases for users to experiment with.
Here are the cases we currently support:

* em\_b\_wave
* em\_convrad
* em\_esmf\_exp
* em\_fire
* em\_grav2d\_x
* em\_heldsuarez
* em\_hill2d\_x
* em\_les
* em\_quarter\_ss
* em\_real
* em\_scm\_xy
* em\_seabreeze2d\_x
* em\_squall2d\_x
* em\_squall2d\_y
* em\_tropical\_cyclone

## Is it possible to run any kind of pre-processing with WRF?

Absolutely! We've included and compiled the WRF Pre-Processing System (WPS),
which you can find in the `/WRF/WPS` directory. The pre-processing utilities are
located in `/WRF/WPS/util`, and since this directory is added to the system
`PATH`, you can run the tools directly from anywhere.

Here’s a list of the available pre-processing tools:

* avg\_tsfc.exe
* calc\_ecmwf\_p.exe
* g1print.exe
* g2print.exe
* height\_ukmo.exe
* int2nc.exe
* mod\_levs.exe
* rd\_intermediate.exe

## How does the process of selecting a use case work?

We’ve compiled all available use cases, and your simulation will run inside the
folder corresponding to the selected case. For example, if you choose the
`em_fire` case (available at [WRF GitHub – em\_fire](https://github.com/wrf-model/WRF/tree/master/test/em_fire))
, your simulation will be executed within that directory.

By default, the simulation will use the example files provided in that folder.
However, if you send us custom input files, we’ll use yours instead.

For instance, the `em_fire` case includes a file named `input_sounding_rain`.
If you provide a file with the same name in your input set, it will override the
default version.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
