# Computational Resources Janitor

Inductiva API has at all times a janitor running in the background that will clean
up any resources that are not being used. The focus of this janitor is to make
sure that no computational resources are being wasted or left idle, for example,
when a user forgets to terminate their machines.

Still, the janitor is mindful of the user's work and will not terminate any
computational resource right away. Instead, it will wait for a certain time
before terminating the resource.

The janitor will also terminate any computational resource that is running for
a certain time, for safeguards.

The janitor follows the following rules before terminating the resources. The
default values are applied when omitted by the user.
- **Total time after launch of computational resource**, independent of any
simulation still running: 7 days
- **Total time of inactivity allowed**, starting from the moment no simulations
are active and resets if any task arrives: 3 minutes

Please note that at the moment, the janitor doesn't preserve the data of running
simulations. So please be mindful of the time your simulations may take to run,
and if this doesn't suffice for your needs, please
[contact us](mailto:support@inductiva.ai).


The allowed machine group maximum validity (`auto_terminate_ts` or
`auto_terminate_minutes`) and the total time of inactivity allowed
(`max_idle_time`) can be defined when initializing the machine group (this is
valid for all machine group types). `max_idle_time` can be a
`datetime.timedelta` object or an integer representing the number of minutes.

```python
import inductiva
import datetime

machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-16",
    data_disk_gb=20,
    max_idle_time=1,
    # or max_idle_time=datetime.timedelta(minutes=1),
    auto_terminate_ts=datetime.datetime.now(datetime.timezone.utc) +
    datetime.timedelta(hours=10),
)
```
