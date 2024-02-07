# User Quotas

Currently, users can launch several computational resources up to a certain limit,
which must satisfy the following quotas:

+ **maximum allowed machines**: 10 machines
+ **maximum allowed vCPUs per resource**: 80 vCPUs
+ **maximum allowed data disk**: 80 GB
+ **maximum allowed price**: 2 USD per hour ($/h)

These quotas are in place over all of the computational resources launched by the user.
That is, you can have two computational resources up at the same time, but, for example, the total number of machines can not exceed 10. Note that for the elastic machine
groups the quotas are applied with the maximum number of machines that can be active at a certain point.

**Maximum allowed machines:** 

Launching a machine group with 11 preemptible machines of the type `c2d-highcpu-4`
contains 44 vCPUs and costs 0.341 $/h. However, it exceeds the allowed limit and
running the following code snippet
```python
mg = inductiva.resources.MachineGroup(
    machine_type="c2d-highcpu-4",
    num_machines=11,
    spot=True
)
```
raises the following error
```bash
Registering MachineGroup configurations:
Registering machine group failed with exception (403)
Reason: Forbidden
HTTP response headers: HTTPHeaderDict({'content-type': 'application/json', 'X-Cloud-Trace-Context': '13cc564cc46afa047fd298749de7f539', 'Date': 'Wed, 07 Feb 2024 14:18:35 GMT', 'Server': 'Google Frontend', 'Content-Length': '96'})
HTTP response body: b'{"detail":"Maximum allowed number of machines is exceeded. Maximum allowed: 10. Requested: 11."}'
```

**Maximum allowed vCPUs:** 

The API provides immediate acess to performant machines with a high number of vCPUs.
However, at the moment, there is a limit to the maximum number of vCPUs that can be
active. This encapsulates two different scenarios:
- launching a machine group with a single machine:
```python
# Price: 1.053 $/h
mg = inductiva.resources.MachineGroup(
    machine_type="c2d-standard-112",
    spot=True
)
```
- launching two machines, one within the limit and another that exceeds it:
```python
# Price: 0.527 $/h
mg = inductiva.resources.MachineGroup(
    machine_type="c2d-standard-56",
    spot=True
)
mg.start()

# Total price: 0.828 $/h
# Total vCPUs: 88
mg2 = inductiva.resources.MachineGroup(
    machine_type="c2d-standard-32",
    spot=True
)
```


This script will fail when launching the second machine group as 10
instances have already been launched.

**Max Allowed Price** refers to the price per hour of all machine
  groups. In this case, a value of 2 USD/h means that resources that
  cost more than 2 USD/h per hour cannot be requested. This
  accumulates over resources, that is, if a machine group that costs 2
  USD/h is already up and running no more resources can be requested.

**Max Allowed Disk Size**: As the name suggests this refers to the
  size of the disks **per machine**. Unlike all the other quotas, this
  does not accumulate over the machines launched.