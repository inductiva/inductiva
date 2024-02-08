# User Quotas

Currently, users can launch several computational resources up to a certain limit,
which must satisfy the following quotas:

+ **maximum allowed price**: 2 USD per hour ($/h)
+ **maximum allowed machines**: 20 machines
+ **maximum allowed vCPUs per resource**: 120 vCPUs
+ **maximum allowed data disk**: 80 GB

These quotas are in place over all of the computational resources launched by the user.
That is, you can have two computational resources up at the same time, but, for
example, the total number of machines can not exceed 10. Note that for the elastic
machine group the quotas are applied with the maximum number of machines that can
be active at a certain point.

**Maximum allowed price:**

The quota limit excess is achieved for example a single machine of type `c2d-standard-56`:
```python
import inductiva

mg = inductiva.resources.MachineGroup("c2d-standard-56")
```
```bash
Error 
```

However, users can take advantage of the computational resource feature that
enable them to select preemptible machines and they may be able to get a
high-performant machine, like a `c2d-highcpu-112`, which fulfils all quotas:
```python
# Cost: 0.870 $/h
mg = inductiva.resources.MachineGroup("c2d-highcpu-112", spot=True)
```

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

When trying to register the second group the following error is raised:
```bash
Registering MachineGroup configurations:
Registering machine group failed with exception (403)
Reason: Forbidden
HTTP response headers: HTTPHeaderDict({'content-type': 'application/json', 'X-Cloud-Trace-Context': 'a7e5789c7404b9bd53ef27f98542b2a5', 'Date': 'Wed, 07 Feb 2024 14:45:17 GMT', 'Server': 'Google Frontend', 'Content-Length': '83'})
HTTP response body: b'{"detail":"Maximum allowed cores is exceeded. Maximum allowed: 80. Requested: 88."}'
```

**Maximum allowed data disk:**

The last quota sets the maximum data disk that users can set for their computational
resources. Here, the quota works per machine, which means you can have a machine group
with four machines all with `100` GB of disk storage dedicated to simulation data.
Only when exceeding this value per machine an error is raised:

```python3
import inductiva

mg = inductiva.resources.MachineGroup("c2-standard-4", data_disk_gb=120)
```


If any of these quotas establish a limit for what you can achieve with Inductiva API,
please [reach out to us](support@inductiva.ai) and we can better understand your
needs.
