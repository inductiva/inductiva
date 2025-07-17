# inductiva **resources** [\[subcommands\]](#subcommands) [\[flags\]](#flags)
The `inductiva resources` command provides utilities for managing computational resources. It allows users to estimate costs, check available machine types, list their active resources, and terminate them.

## Subcommands
### `start`
Start computational resources.

```bash
inductiva resources start
```

### `cost` [\[flags\]](#flags-for-cost)
Estimate the hourly cost of a machine or group of machines in the cloud.

```bash
inductiva resources cost <MACHINE_TYPE>
```

Estimate the cost of using 4 spot machines of type **c2-standard-8**:

```bash
$ inductiva resources cost c2-standard-8 --spot -n 4
Estimated total cost (per machine): 0.445 (0.111) $/h.
```

<h4 id="flags-for-available">Flags</h4>

**`--spot`** (default: `false`)

Estimate the cost for a **spot instance** of the specified machine type.

---

**`--num-machines, -n`** (default: 1)

Number of machines to include in the cost estimate.

---

**`--zone`** (default: `europe-west1-b`)

Zone where the machine(s) will be launched. Pricing may vary by zone.

### `info`
Print information related to a machine group.

```bash
inductiva resources info <MACHINE_GROUP_NAME>
```

### `available` [\[flags\]](#flags-for-available)
List available machine types.

```bash
inductiva resources available
```

List the machines of the `c3d` family.

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

<h4 id="flags-for-available">Flags</h4>

**`--family, -f`**
Filter available machine types by CPU family (e.g., `c3d`, etc.).
Supports one or more values.

---

**`--provider, -p`** (default: `GCP`)

Filter the available types by provider. Options:
- `LOCAL`: list locally available machine types.
- `GCP`: list machine types available in GCP

### `list ls`
List all active computational resources associated with your account.

This subcommand provides a complete snapshot of your currently running resources, including details such as machine type, whether it's elastic or spot-based, the number of active machines, and pricing.

```bash
inductiva resources list
```

or using the alias:

```bash
inductiva resources ls
```

> Note: No resources will be shown if none are currently active.

Sample output:

```sh
$ inductiva resouces ls
  Active Resources:

 NAME                            MACHINE TYPE     ELASTIC     TYPE       # MACHINES     DATA SIZE IN GB     SPOT     CREATED AT (UTC)     IDLE TIME      MAX COST ($/HOUR)
 api-sjv5giwarf49kt89f5em1z3rv   c4-highcpu-4     False       standard   0/1            10                  True     17/07, 20:59:32      None/0:03:00   0.689884
```

### `terminate` [\[flags\]](#flags-for-terminate)
Terminate one or more active computational resources.

```bash
inductiva resources terminate <MACHINE_GROUP_NAME> [<MACHINE_GROUP_NAME> ...]
```

<h4 id="flags-for-available">Flags</h4>

**`--yes, y`**

Automatically answers **yes** to all confirmation prompts, allowing operations to proceed without interactive confirmation. Use with caution for irreversible actions.

---

**--all, -a**

Terminates all resouces.

### Estimating Machine Cost

You can estimate the costs of the computational resources you 
plan to use per hour. The CLI provides a cost estimation tool
that considers the machine type, usage duration, and number of
machines you wish to use.

## Flags
### `-h, --help`

Show help message and exit.

## Need Help?
Run the following command for more details:

```sh
inductiva resources --help
```

