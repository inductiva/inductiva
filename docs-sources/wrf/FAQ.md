**Find answers to commonly asked questions about WRF.**

<br>

# FAQ

## In your tutorials, you use the case `em_real`. What other cases are available?
WRF includes a variety of pre-configured test cases for users to explore. 
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

<br>

## Is it possible to run any kind of pre-processing with WRF?
Absolutely! We include and compile the WRF Pre-Processing System (WPS), 
located in the `/WRF/WPS` directory. The pre-processing utilities reside in `/WRF/WPS/util`, 
and since this directory is added to the system `PATH`, you can run these tools directly from anywhere.

Here’s a list of the available pre-processing tools:

* avg\_tsfc.exe
* calc\_ecmwf\_p.exe
* g1print.exe
* g2print.exe
* height\_ukmo.exe
* int2nc.exe
* mod\_levs.exe
* rd\_intermediate.exe

<br>

## How does the process of selecting a use case work?
All available use cases are compiled, and your simulation runs inside the folder corresponding 
to the chosen case. For example, selecting the `em_fire` case (available at [WRF GitHub – em\_fire](https://github.com/wrf-model/WRF/tree/master/test/em_fire)) means your simulation will execute within that directory.

By default, the simulation uses the example files provided in the folder. However, if you provide custom input files, those will override the defaults.

For instance, the `em_fir`e` case includes a file named `input_sounding_rain`. If you provide a file with the same name in your input set, it will replace the default version.

<br>

## Why is my `gen_gif.py` command failing?
If you encounter an error like this:

```
urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>

ERROR conda.cli.main_run:execute(125): `conda run python /scripts/gen_gif.py --files ... --output-dir . --fps 3 --var RAINNC` failed. (See above for error)
INFO:    Cleanup error: while stopping driver for /var/lib/apptainer/mnt/session/final: fuse-overlayfs exited: fuse: reading device: Software caused connection abort
```

It likely means that **Cartopy is trying to download map data**, but your
machine **does not have an internet connection**. Without access to these
external resources, the script fails with the error shown above.

*Please inform us if this occurs*. We can update the simulator to include the necessary Cartopy files 
so they are cached locally and no internet access is needed during the GIF generation process.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
