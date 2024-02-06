# User Resource Quotas

The `inductiva` package comes with the following quotas for
computational resources:

+ **Max Allowed Instances**: 10 machines
+ **Max Allowed Price**: 2 USD
+ **Max Allowed Cores**: 80 cores
+ **Max Allowed Disk Size**: 80 GB

Where **Max Allowed Instances** controls the number of instances that
a user can have active at any point. For example, the next code
snippets will both raise an error:

```python
mg = inductiva.resources.MachineGroup(
    machine_type="c2-standard-8",
    num_machines=11,
)
mg.start()
```

This will fail launching the machine group since more machines than
those allowed are being requested.

```python
mg = inductiva.resources.MachineGroup(
    machine_type="c2-standard-8",
    num_machines=10,
)
mg.start()

mg = inductiva.resources.MachineGroup(
    machine_type="c2-standard-8",
    num_machines=1,
)
mg.start()
```

This script will fail when launching the second machine group as 10
instances have already been launched.

**Max Allowed Price** refers to the price per hour of all machine
  groups. In this case, a value of 2 USD means that resources that
  cost more than 2 USD per hour cannot be requested. This accumulates
  over resources, that is, if a machine group that costs 2 USD is
  already up and running no more resources can be requested.

**Max Allowed cores** works similarly to **Max Allowed
  Instances**. That is, the total number of cores in all active
  machines groups cannot exceed 80 cores.

**Max Allowed Disk Size**: As the name suggests this refers to the
  size of the disks **per machine**. Unlike all the other quotas, this
  does not accumulate over the machines launched.