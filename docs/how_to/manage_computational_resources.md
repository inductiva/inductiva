# Manage Computational Resources

Once you have launched your computational resources, there are a few API methods
that help manage them and see their status. Let's go over them one by one.

## Get your active computational resources

In case computational resources have been launched in another Python session and
want to manage or re-use them in another script, users can fetch the respective
instance via the `get` method as follows:

```python
## Obtain a list with instances of all active computational resources
>>> resources_list = inductiva.resources.machine_groups.get()
>>> print(resources_list)
[MPICluster(name="api-23zssj6oq77xxsot3o0nhax3d"),
 ElasticMachineGroup(name="api-45fetsr58okcs0x6j9m0vsi2z"),
 MachineGroup(name="api-4kken08fnoxuu5zjjak6ak2xe")]
```

## List your active computational resources

When you just want to check the active resources you can quickly list the
information of each one either via Python or via the CLI.

**Python**
```python
inductiva.resources.machine_groups.get()
```

**CLI**
```bash
$ inductiva resources list
```

One obtains for example the following information:
```
Active Resources:

 NAME                            MACHINE TYPE     ELASTIC     TYPE       # MACHINES     DATA SIZE IN GB     SPOT     STARTED AT (UTC)     IDLE TIME         MAX COST ($/HOUR)
 api-3ejvh64mxuxnfcv3yxdhoyjuj   c2-standard-4    False       standard   5/5            50                  False    10 Jul, 16:23:00     0:04:15/0:30:00   1.4909
 api-5014txg0rwx3jbbpf6y0ndzmv   c2d-highmem-16   False       mpi        3/3            10                  False    10 Jul, 16:22:04     0:05:12/0:30:00   3.27774
 api-es9sjockjymvkwfmjioibfw8p   c2-standard-8    False       standard   2/2            60                  False    10 Jul, 16:23:25     0:03:50/0:30:00   1.08312
```

## Terminate the active computational resources

When you have finished using your computational resources, don't forget to terminate
them. An advantage of Inductiva API is to control your computational
resources and avoid them being idle.

Hence, you can either terminate your computational resources via Python or the CLI. 

In Python, you will need to instantiate the computational resource object, for
example, with the `get` method and then do `machine.terminate()` for example. 

Via CLI, this process is much easier and you can either terminate a machine on
at a time:
```bash
inductiva resources terminate api-45fetsr58okcs0x6j9m0vsi2z
```

Or terminate them all at once. You will need to confirm this action:
```bash
inductiva resources terminate --all
```

Both are a blocking call that will only finish when the machines have terminated,
in this way no computational resources are left up.
