# Resources

As you might know by now, Inductiva API provides a simple way to [launch
dedicated resources]() where you can run your simulations. Before launching any
resources, users can use the CLI to gather more information about the right
resources to launch with the `available` and `cost` subcommands.

With them, user can list all the available machine types together with details,
```bash
$ inductiva resources available
Available machine types

machine-type: [cores-available]
c2-standard-: [4, 8, 16, 30, 60]
c3-standard-: [4, 8, 22, 44, 88, 176]
c2d-standard-: [2, 4, 8, 16, 32, 56, 112]
c2d-highcpu-: [2, 4, 8, 16, 32, 56, 112]
e2-standard-: [2, 4, 8, 16, 32]
n2-standard-: [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]
n2d-standard-: [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]
n1-standard-: [1, 2, 4, 8, 16, 32, 64, 96]
e2-highcpu-: [2, 4, 8, 16, 32]

 E.g. of machine-type: c2-standard-8
```

and thereafter, use one of the machine types to get an estimate of the cost
per hour to use it:
```bash
$ inductiva resources cost c2-standard-8 --spot -n 4
Estimated total cost (per machine): 0.445 (0.111) $/h.
```

When you have already decided and launch a few computational resources, you can
use the list subcommand to get an overview of the resources you have running:
```bash
$ inductiva resources list
Active Resources:
Name                           Machine Type    Elastic    Type        # machines    Disk Size in GB  Spot    Started at (UTC)
-----------------------------  --------------  ---------  --------  ------------  -----------------  ------  ------------------
api-abs7c8hwccu0v6cppo9nc40z6  c2-standard-30  False      standard             5                 70  True    01 Feb, 23:06:52
api-mym74u9i9xppi9ngzze3xueae  c2-standard-8   False           mpi             4                 70  True    02 Feb, 00:29:53
```

Finally, the CLI also allows to terminate the resources that are no longer required with
the `terminate` subcommand. Users can either choose a specific resource with the
`--name` flag or terminate all the resources with the `--all` flag. Any of the steps
require confirmation from the user before proceeding.
```bash
$ inductiva resources terminate --all
Confirm the termination of all active machines? (y/[N]) y
Terminating all active computational resources.
Terminating MachineGroup(name="api-abs7c8hwccu0v6cppo9nc40z6"). This may take a few minutes.
Machine Group api-abs7c8hwccu0v6cppo9nc40z6 with c2-standard-30 machines successfully terminated in 0:01:07.
Terminating MachineGroup(name="api-mym74u9i9xppi9ngzze3xueae"). This may take a few minutes.
Machine Group api-mym74u9i9xppi9ngzze3xueae with c2-standard-8 machines successfully terminated in 0:01:12.
```
