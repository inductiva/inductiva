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

Moreover, to avoid cluttering the shared pool of resources, which is mostly meant for
testing purposes, the janitor will terminate any task that is running also for a
certain time. 

The janitor follows the following rules before terminating the resources:
- **Total time after launch of computational resource**, independent of any
simulation still running: 36h 
- **Total time of inactivity allowed**, starting from the moment no simulations are
on queue and reset if any arrive: 30 min 
- **Total time of task on default machine group**: 4h

Please note that at the moment, the janitor is clumsy and doesn't know how to handle
the data of running simulations. So please be mindful of the time your simulations
may take to run, and if this doesn't suffice for your needs, please [contact us](support@inductiva.ai).

