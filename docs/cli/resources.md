# Resources

As you might know by now, Inductiva API provides a simple way to [launch dedicated resources](../how_to/computational_resources.md)
where you can run your simulations. Before launching any
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

       NAME                                MACHINE TYPE         ELASTIC         TYPE           # MACHINES         DATA SIZE IN GB         SPOT         STARTED AT (UTC)
       api-p3kun5wyta1hacstu4xk38ujr       c2-standard-8        False           mpi            2                  10                      False        08 Feb, 12:59:10
       api-rdqprn82417bsd7id1qnac4c6       c2-standard-4        False           standard       16                 10                      False        08 Feb, 12:58:28
```

Finally, the CLI also allows the termination of resources that are no longer required with
the `terminate` subcommand. Users can either choose a specific resource by
providing its name or terminate all the resources with the `--all` flag. Any of the steps
require confirmation from the user before proceeding. Here, we choose to terminate
all the resources:

```bash
$ inductiva resources terminate --all
You are about to terminate ALL resources.
Are you sure you want to proceed (y/[N])? y
Terminating MPICluster(name="api-p3kun5wyta1hacstu4xk38ujr"). This may take a few minutes.
MPI Cluster api-p3kun5wyta1hacstu4xk38ujr with c2-standard-8 x2 machines successfully terminated in 0:01:10.
Terminating MachineGroup(name="api-rdqprn82417bsd7id1qnac4c6"). This may take a few minutes.
Machine Group api-rdqprn82417bsd7id1qnac4c6 with c2-standard-4 machines successfully terminated in 0:01:18.
```
