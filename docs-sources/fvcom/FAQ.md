**Find answers to commonly asked questions about FVCOM.**

<br>

# FAQ

## 1. Where can I find the configuration settings used to compile each binary?
The configuration settings used in each binary can be found in the respective `make.inc` files:
- For **fvcom**, refer to the [make.inc](https://github.com/inductiva/kutu/blob/main/simulators/fvcom/v5.1.0/make.inc).
- For **fvcom_estuary**, refer to the [make_estuary.in](https://github.com/inductiva/kutu/blob/main/simulators/fvcom/v5.1.0/make_estuary.inc).

<br>

## 2. What should I do if I encounter issues with the input namelist file?
If you encounter issues with the input namelist file, use the `create_namelist` parameter to generate a valid `.nml` file in your working directory. To do this, modify the following section of the FVCOM code as shown below:

```python
# Run simulation 
task = fvcom.run(
    input_dir=input_dir,
    working_dir="run/",
    create_namelist="tst",
    n_vcpus=1,
    on=machine_group)
```

You can now download your outputs and use the newly created `.nml` file as a
base of your simulation and start working from this valid file.

<br>

## 3. How can I debug a FVCOM simulation?
To obtain more detailed information in the `stdout` and `stderr` output files, you can enable the `debug` parameter. To do this, modify the following section of the FVCOM code as follows:

```python
# Run simulation
task = fvcom.run(
    input_dir=input_dir,
    working_dir="run/",
    case_name="tst",
    n_vcpus=1,
    debug=7,
    on=machine_group)
```

The `debug` parameter can be set to various levels, each providing different levels of detail:
debg=0 → Default debug log 
debg=1 → Debug IO filenames 
debg=2 → Debug scalar values 
debg=4 → Debug subroutine names 
debg=5 → Debug subroutine IO 
debg=6 → Debug vector values 
debg=7 → Full debug information 

<br>

## 4. What should I do if I encounter issues with the timezone argument in the namelist file?
If you encounter issues with the timezone argument in the `.nml` file, try setting it to `None` or `UTC` as a workaround. For more details, refer to the [bug report](https://github.com/FVCOM-GitHub/FVCOM/issues/27).

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
