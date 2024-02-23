# Manage your Computational Resources

The Inductiva API provides a simple way to [launch dedicated resources](../how_to/computational_resources.md) for running your simulations. With the CLI, you 
can effectively manage these computational resources, from selection and launch to 
termination, directly from your terminal. For further details on each command and 
additional options, refer to the CLI's built-in help system using the `--help` flag.

## Discover Available Resources

Before launching any resources, you can use the CLI to go through the variety 
of machine types available and their associated costs using the `available` and 
`cost` subcommands:

```bash
$ inductiva resources available
Machine types provided in Google Cloud

# You would get an output listing all available machines:

c2: Intel Xeon Cascade Lake (2nd Gen) processor.
  > c2-standard-  [2, 4, 8, 16, 32, 60]                         

c3: Intel Xeon Sapphire Rapids (4th Gen) processor.
  > c3-highcpu-   [4, 8, 22, 44, 88, 176]                       
  > c3-standard-  [4, 8, 22, 44, 88, 176]                       
  > c3-highmem-   [4, 8, 22, 44, 88, 176]                       

h3: (Available Soon) Intel Xeon Sapphire Rapids (4th Gen) processor.
Simultaneous multithreading disabled, i.e., vCPU represents an entire core.
  > h3-standard-  [88]                                          

c2d: AMD EPYC Milan (3rd Gen) processor.
  > c2d-highcpu-  [2, 4, 8, 16, 32, 56, 112]                    
  > c2d-standard- [2, 4, 8, 16, 32, 56, 112]                    
  > c2d-highmem-  [2, 4, 8, 16, 32, 56, 112]                    

c3d: AMD EPYC Genoa (4th Gen) processor.
  > c3d-highcpu-  [4, 8, 16, 30, 60, 90, 180, 360]              
  > c3d-standard- [4, 8, 16, 30, 60, 90, 180, 360]              
  > c3d-highmem-  [4, 8, 16, 30, 60, 90, 180, 360]              

e2: Intel Xeon (up to Skylake, 1st Gen) and AMD EPYC (up to Milan, 3rd Gen) processors.
Automatically selected based on availability.
  > e2-highcpu-   [2, 4, 8, 16, 32]                             
  > e2-standard-  [2, 4, 8, 16, 32]                             
  > e2-highmem-   [2, 4, 8, 16]                                 

n2: Intel Xeon Ice Lake and Cascade Lake processors (3rd and 2nd Gen).
Cascade Lake default up to 80 vCPUs and Ice Lake for larger machines.
  > n2-highcpu-   [2, 4, 8, 16, 32, 48, 64, 80, 96]             
  > n2-standard-  [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]        
  > n2-highmem-   [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]        

n2d: AMD EPYC Milan or ROME processors (3rd and 2nd Gen).
  > n2d-highcpu-  [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]   
  > n2d-standard- [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]   
  > n2d-highmem-  [2, 4, 8, 16, 32, 48, 64, 80, 96]             

n1: Intel Xeon (up to Skylake, 1st Gen) processor.
Automatically selected based on availability.
  > n1-highcpu-   [1, 2, 4, 8, 16, 32, 64, 96]                  
  > n1-standard-  [1, 2, 4, 8, 16, 32, 64, 96]                  
  > n1-highmem-   [1, 2, 4, 8, 16, 32, 64, 96]  
```

If you want more detailed information you can use the `-v` flag, and you can focus
on a specific series by using the `-s` flag. 

For example:

```bash
$ inductiva resources available -s c3d -v
Machine types provided in Google Cloud
c3d: AMD EPYC Genoa (4th Gen) processor.
  > c3d-highcpu-  [4, 8, 16, 30, 60, 90, 180, 360]              -> 2 GB of memory per vCPU.
  > c3d-standard- [4, 8, 16, 30, 60, 90, 180, 360]              -> 4 GB of memory per vCPU and possible local ssd integration.
  > c3d-highmem-  [4, 8, 16, 30, 60, 90, 180, 360]              -> 4 GB of memory per vCPU.
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