# User Quotas

In the current version of the API (0.7) users can launch several computational resources
up to a certain limit, which must satisfy the quotas we present on this section.
The quotas are in place to establish a limit over the totality of the computational
resources launched by the user. 

Note that for the Elastic Machine Groups the quotas are applied to the maximum
number of machines requested, not the number of machines that are currently
active.

## Maximum allowed price: 2 $/h

The quota limit can be reached even with a single machine. For example, if you
select a machine of type `c2d-standard-56`, which has a cost of 2.799 $/h, this
is what will happen:

```python
import inductiva

mg = inductiva.resources.MachineGroup("c2d-standard-56")
mg.start()
```
```bash
Registering MachineGroup configurations:
> Name:         api-ikg3n235fcxnkh7u97etiv93n
> Machine Type: c2d-standard-56
> Data disk size:    10 GB
> Number of machines: 1
> Spot:               False
> Estimated cloud cost of machine group: 2.799 $/h
Starting MachineGroup(name="api-bbncoqy189emgs1768qxbucaf"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Starting machine group failed with exception (400)
Reason: Bad Request
HTTP response headers: HTTPHeaderDict({'content-type': 'application/json', 'X-Cloud-Trace-Context': 'e38080d4b070af7fb9efab64f25f9c87', 'Date': 'Fri, 09 Feb 2024 13:28:54 GMT', 'Server': 'Google Frontend', 'Content-Length': '90'})
HTTP response body: b'{"detail":"Quota exceeded: cost ($/h)\\nRequested: 2.79927\\nIn use: 0.0\\nMax allowed: 2\\n"}'
```

However, users can take advantage of **spot instances** that, despite being
preemptible by the cloud provider, are much less expensive. For example, you
can get an even better machine than the one above, if you select it to run as
spot instance:
```python
# Cost: 0.870 $/h
import inductiva

# Launching a c2d-highcpu-112 machine as spot instance
mg = inductiva.resources.MachineGroup("c2d-highcpu-112", spot=True)
mg.start()
```
```bash
Registering MachineGroup configurations:
> Name:         api-io3bspoyy13astzyih291c2m8
> Machine Type: c2d-highcpu-112
> Data disk size:    10 GB
> Number of machines: 1
> Spot:               True
> Estimated cloud cost of machine group: 0.870 $/h
Starting MachineGroup(name="api-ibg2s91an4po453ry7cwv8ie9"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Machine Group api-ibg2s91an4po453ry7cwv8ie9 with c2d-highcpu-112 machines successfully started in 0:00:20.
```

Notice, that you should also take into account the machine **variant** you
choose. For example, in the case of the CPU series `c2d`, machines of the 
`highmem` variant (8 GB RAM per vCPU) are more expensive than machines of
the `standard` variant (4 GB RAM per vCPU), which, in turn, are more expensive
than machines of the `highcpu` variant (2 GB RAM per vCPU). The amount of RAM
attached to a machine is one of the most determining aspects of its price.

## Maximum allowed machines: 20

Launching a machine group with 20 preemptible machines of the type 
`c2d-highcpu-4` contains 44 vCPUs and costs in total 0.341 $/h. However, it
exceeds the allowed limit of machines and running the following code snippet
raises the following error:
```python
import inductiva

mg = inductiva.resources.MachineGroup(
    machine_type="c2d-highcpu-4",
    num_machines=21,
    spot=True
)

mg.start()
```
raises the following error
```bash
Registering MachineGroup configurations:
> Name:         api-aizo46wdj6eiiqostiidxwha1
> Machine Type: c2d-highcpu-4
> Data disk size:    10 GB
> Number of machines: 21
> Spot:               True
> Estimated cloud cost of machine group: 0.652 $/h
Starting MachineGroup(name="api-aizo46wdj6eiiqostiidxwha1"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Starting machine group failed with exception (400)
Reason: Bad Request
HTTP response headers: HTTPHeaderDict({'content-type': 'application/json', 'X-Cloud-Trace-Context': '5f4cec0da72dc3bbfda15c43aed4393f', 'Date': 'Fri, 09 Feb 2024 12:51:51 GMT', 'Server': 'Google Frontend', 'Content-Length': '92'})
HTTP response body: b'{"detail":"Quota exceeded: number of machines\\nRequested: 21\\nIn use: 0\\nMax allowed: 20\\n"}'
```

## Maximum allowed vCPUs: 240

The quota limit of the maximum allowed vCPUs is harder to reach with a single
machine group. In most cases, the two above quotas will be reached before this
one.

Still, there is an example where the limit can be reached and uses the ability
of having several machine groups active at the same time.
```python
import inductiva

# Sucessfull launch with cost 1.045 $/h and total 128 vCPUs
mg1 = inductiva.resources.MachineGroup(
    "e2-highcpu-32",
    num_machines=4,
    spot=True)
mg1.start()

# Sucessfull launch. Total cost 1.959 $/h and total 240 vCPUs 
mg2 = inductiva.resources.MachineGroup(
    "e2-highcpu-16",
    num_machines=7,
    spot=True)
mg2.start()

# Failed launch. Total cost 1.975 $/h, total machines 11 and total vCPUs 242.
mg3 = inductiva.resources.MachineGroup(
    "e2-highcpu-2", 
    num_machines=1,
    spot=True)
mg3.start()
```
The following error is raised on the last start:
```bash
Starting MachineGroup(name="api-cs82oayc9rxgmig8zzmg4mh59"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Starting machine group failed with exception (400)
Reason: Bad Request
HTTP response headers: HTTPHeaderDict({'content-type': 'application/json', 'X-Cloud-Trace-Context': '05b4aeec3827418341021ade1e4e54a3', 'Date': 'Fri, 09 Feb 2024 13:10:03 GMT', 'Server': 'Google Frontend', 'Content-Length': '91'})
HTTP response body: b'{"detail":"Quota exceeded: number of cores\\nRequested: 2\\nIn use: 240\\nMax allowed: 240\\n"}'
```

## Maximum allowed data disk: 20 GB

The last quota sets the maximum data disk that users can set for their 
computational resources. Here, the quota works per machine. The current limit
on disk space is 20 GB, so the next request will hit the limit:

```python
import inductiva

mg = inductiva.resources.MachineGroup("c2-standard-4", data_disk_gb=30)
```
```bash
Registering MachineGroup configurations:
> Name:         api-0w5epve63qbirww7j80f4suqw
> Machine Type: c2-standard-4
> Data disk size:    30 GB
> Number of machines: 1
> Spot:               False
> Estimated cloud cost of machine group: 0.230 $/h
Starting MachineGroup(name="api-0w5epve63qbirww7j80f4suqw"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Starting machine group failed with exception (400)
Reason: Bad Request
HTTP response headers: HTTPHeaderDict({'content-type': 'application/json', 'X-Cloud-Trace-Context': '39b0e0ed652420807dba1c2fd471801c', 'Date': 'Fri, 09 Feb 2024 12:48:37 GMT', 'Server': 'Google Frontend', 'Content-Length': '76'})
HTTP response body: b'{"detail":"Max disk size exceeded.\\nRequested: 30 GB\\nMax allowed: 20 GB\\n"}'
```

**Maximum allowed vCPUs:** 

The API provides immediate acess to performant machines with a high number of
vCPUs. However, at the moment, there is a limit to the maximum number of vCPUs
that can be active. This encapsulates two different scenarios:

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

The last quota sets the maximum data disk that users can set for their
computational resources. Here, the quota works per machine, which means you can
have a machine group with four machines all with `100` GB of disk storage
dedicated to simulation data. Only when exceeding this value per machine an
error is raised:

```python
import inductiva

mg = inductiva.resources.MachineGroup("c2-standard-4", data_disk_gb=120)
```


If any of these quotas establish a limit for what you can achieve with Inductiva
API, please [reach out to us](mailto:support@inductiva.ai) and we can better
understand your needs.

## How to monitor quotas

Inductiva's cli provides an easy way to monitor your quotas usage. Simply use 
the command `inductiva quotas list`. This will output a detailed list of all
the quotas together with their current usage. For example:

```bash
$ inductiva quotas list

       NAME                               IN_USE         MAX_ALLOWED
       total_num_vcpus                    0              240
       total_num_machines                 0              20
       cost_per_hour                      0              2
       machine_disk_size_gb               n/a            20
       machine_group_idle_minutes         n/a            30
       machine_group_lifetime_hours       n/a            36
```

Alternatively, you can also get the quotas directly from our python client:

```python
import inductiva

inductiva.users.get_quotas()
```
