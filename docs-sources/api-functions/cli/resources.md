# resources

## Overview
The `inductiva resources` command provides utilities for managing computational resources. It allows users to estimate costs, check available machine types, list their active resources, and terminate them.

## Usage
```bash
inductiva resources [-h] {cost,info,available,list} ...
```

## Options
- `-h, --help`  
  Show this help message and exit.

## Available Subcommands
### `start`
Start computational resources.

#### Usage
```bash
inductiva resources start
```

### `cost`
Estimate the cost of a machine in the cloud.

#### Usage
```bash
inductiva resources cost
```

### `info`
Print information related to a machine group.

#### Usage
```bash
inductiva resources info
```

### `available`
List available machine types.

#### Usage
```bash
inductiva resources available
```

### `list` (alias: `ls`)
List currently active resources.

#### Usage
```bash
inductiva resources list
```

or using the alias:

```bash
inductiva resources ls
```

### `terminate`
Terminate resources.

#### Usage
```bash
inductiva resources terminate
```

## Examples

### Listing Available Machine Types
```bash
$ inductiva resources available

MACHINE TYPE            VCPUS     GPUS                     MEMORY (GB)     PRICE/HOUR (USD)     ZONE
c2-standard-4           4         n/a                      16              0.2297               europe-west1-b
c2-standard-8           8         n/a                      32              0.4594               europe-west1-b
c2-standard-16          16        n/a                      64              0.9188               europe-west1-b
c2-standard-30          30        n/a                      120             1.72275              europe-west1-b
c2-standard-60          60        n/a                      240             3.4455               europe-west1-b
c2d-highcpu-2           2         n/a                      4               0.082538             europe-west1-b
c2d-highcpu-4           4         n/a                      8               0.165076             europe-west1-b
c2d-highcpu-8           8         n/a                      16              0.330152             europe-west1-b
...
```

### Listing a specific machine Family
This is how you list the machines of the c3d family.

```bash
$ inductiva resources available -f c3d

MACHINE TYPE            VCPUS     GPUS     MEMORY (GB)     PRICE/HOUR (USD)     ZONE
c3d-highcpu-4           4         n/a      8               0.16508132           europe-west1-b
c3d-highcpu-8           8         n/a      16              0.33016264           europe-west1-b
c3d-highcpu-16          16        n/a      32              0.66032528           europe-west1-b
c3d-highcpu-30          30        n/a      59              1.233750645          europe-west1-b
c3d-highcpu-60          60        n/a      118             2.46750129           europe-west1-b
c3d-highcpu-90          90        n/a      177             3.701251935          europe-west1-b
c3d-highcpu-180         180       n/a      354             7.40250387           europe-west1-b
c3d-highcpu-360         360       n/a      708             14.80500774          europe-west1-b
...
```

### Estimating Machine Cost

You can estimate the costs of the computational resources you 
plan to use per hour. The CLI provides a cost estimation tool
that considers the machine type, usage duration, and number of
machines you wish to use.

Consider the following example, where you wish to estimate
the cost of using f4 machines of type **c2-standard-8**:

```bash
$ inductiva resources cost c2-standard-8 --spot -n 4
Estimated total cost (per machine): 0.445 (0.111) $/h.
```

### List Your Active Resources

You can use the `list` subcommand to get an overview of your 
active computational resources:

```bash
$ inductiva resources list
Active Resources:

       NAME                                MACHINE TYPE         ELASTIC         TYPE           # MACHINES         DATA SIZE IN GB         SPOT         STARTED AT (UTC)
       api-p3kun5wyta1hacstu4xk38ujr       c2-standard-8        False           mpi            2                  10                      False        08 Feb, 12:59:10
       api-rdqprn82417bsd7id1qnac4c6       c2-standard-4        False           standard       16                 10                      False        08 Feb, 12:58:28
```

### Terminate Resources

Finally, you can terminate computational resources that are no longer needed
using `terminate` subcommand. You can either choose a specific resource 
by providing its name or terminate all the resources with the `--all` flag.
All tasks running in the resources you are terminating will be immediately killed,
no matter what stage they are in. 
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

## Additional Resources
For more details, visit the [Inductiva API documentation](#).

