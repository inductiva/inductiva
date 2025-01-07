# Manage your Computational Resources

The Inductiva API provides a simple way to [launch dedicated resources](https://tutorials.inductiva.ai/how_to/manage_computational_resources.html) 
for running your simulations. With the CLI, you can effectively manage these 
computational resources, from selection and launch to termination, directly from 
your terminal. For further details on each command and additional options, refer 
to the CLI's built-in help system using the `--help` flag.

## Discover Available Resources

Before launching any resources, you can use the CLI to go through the variety 
of machine types available and their associated costs using the `available` and 
`cost` subcommands:

```bash
$ inductiva resources available
Machine types provided in Google Cloud

# You would get an output listing all available machines:

CPU family: c2
   > c2-standard- [4, 8, 16, 30, 60]                       

CPU family: c3
   > c3-highcpu- [4, 8, 22, 44, 88, 176] 
   > c3-highmem- [4, 8, 22, 44, 88, 176] 
   > c3-standard- [4, 8, 22, 44, 88, 176] (-lssd)                   

...

```

You can focus on a specific series by using the `-f` flag. 

For example:

```bash
$ inductiva resources available -f c3d
Machine types provided in Google Cloud
c3d: AMD EPYC Genoa (4th Gen) processor.
  > c3d-highcpu-  [4, 8, 16, 30, 60, 90, 180, 360]
  > c3d-standard- [4, 8, 16, 30, 60, 90, 180, 360]
  > c3d-highmem-  [4, 8, 16, 30, 60, 90, 180, 360]
```
## Estimate Costs

You can estimate the costs of the computational resources you plan to use per hour. 
The CLI provides a cost estimation tool that considers the machine type, usage duration, 
and number of machines.

Consider the following example, where you wish to estimate the cost of four **c2-standard-8** machines:

```bash
$ inductiva resources cost c2-standard-8 --spot -n 4
Estimated total cost (per machine): 0.445 (0.111) $/h.
```

## List Active Resources

Once you've decided and launched your resources, you can use the `list` subcommand 
to get an overview of your active computational resources:

```bash
$ inductiva resources list
Active Resources:

       NAME                                MACHINE TYPE         ELASTIC         TYPE           # MACHINES         DATA SIZE IN GB         SPOT         STARTED AT (UTC)
       api-p3kun5wyta1hacstu4xk38ujr       c2-standard-8        False           mpi            2                  10                      False        08 Feb, 12:59:10
       api-rdqprn82417bsd7id1qnac4c6       c2-standard-4        False           standard       16                 10                      False        08 Feb, 12:58:28
```

## Terminate Resources

Finally, you can terminate computational resources that are no longer needed through 
the CLI with the `terminate` subcommand. You can either choose a specific resource 
by providing its name or terminate all the resources with the `--all` flag. 
**Any of the steps require user confirmation before proceeding.** 

For example, you can choose to terminate all the resources:

```bash
$ inductiva resources terminate --all
# The CLI always prompts for confirmation before terminating resources
You are about to terminate ALL resources.
Are you sure you want to proceed (y/[N])? y
Terminating MPICluster(name="api-p3kun5wyta1hacstu4xk38ujr"). This may take a few minutes.
MPI Cluster api-p3kun5wyta1hacstu4xk38ujr with c2-standard-8 x2 machines successfully terminated in 0:01:10.
Terminating MachineGroup(name="api-rdqprn82417bsd7id1qnac4c6"). This may take a few minutes.
Machine Group api-rdqprn82417bsd7id1qnac4c6 with c2-standard-4 machines successfully terminated in 0:01:18.
```