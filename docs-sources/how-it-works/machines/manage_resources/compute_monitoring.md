# Automated Resource Monitoring

Inductiva API has at all times a monitoring service running in the background
that will clean up any resources that are not being used. The focus of this
service is to make sure that no computational resources are being wasted or
left idle, for example, when a user forgets to terminate their machines.

Still, this service is mindful of the user's work and will not terminate any
computational resource right away. Instead, it will wait for a certain time
before terminating the resource.

The following rules are verified before terminating the resources, which can be specified by the user:
- **Total time of inactivity allowed**, starting from the moment no simulations
are active and resets if any task arrives: `max_idle_time`. The default value
is 3 minutes.
- **Total time after launch of computational resource**, independent of any
simulation still running: `auto_terminate_ts` or `auto_terminate_minutes`. There
is no default value. If not specified, the machine group will never terminate
while having active tasks assigned.

Please note that at the moment, the background service doesn't preserve the
data of running simulations. So please be mindful of the time your simulations
may take to run when specifying `auto_terminate_ts` or `auto_terminate_minutes`.

The above mentioned attributes can be defined when initializing the resource.

```python
import inductiva
import datetime

machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-16",
    data_disk_gb=20,
    max_idle_time=1,  # terminate after 1 minute of inactitivy
    # or
    # max_idle_time=datetime.timedelta(seconds=60),
    auto_terminate_ts=datetime.datetime.now(datetime.timezone.utc) +
    datetime.timedelta(hours=1),  # terminate 60 minutes after launch
    # or
    # auto_terminate_minutes=60,
)
```
