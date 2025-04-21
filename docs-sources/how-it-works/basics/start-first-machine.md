# Start your first cloud machine with Inductiva

In this guide, we'll take you step-by-step while you start your first cloud machine with Inductiva.
Here's what you'll need to know:

1️⃣ Pick the right cloud machine for your simulation use-case
Learn how to navigate performance and price options among hundreds of cloud machines.

2️⃣ Start your first dedicated cloud machine
Start a machine with 3 lines of Python code without worrying about any configuration details.

3️⃣ Terminate your cloud machine

## 1️⃣ Pick your cloud machine

Inductiva provides you access to hundreds of cutting-edge cloud machines and empowers you to select the best option for your simulation use-case.
Follow <a href="pick_cloud_machine.html">the guide</a> we prepared.

After narrowing down the machine that best fits your needs, let’s move on and get it started!

## 2️⃣ Start your first dedicated cloud machine

Starting a machine with Inductiva means you’re allocating it exclusively for your own use - no queue, and no waiting time to run your tasks.

The following python script allows you to start a cloud machine with Inductiva:

```python
import inductiva
cloud_machine = inductiva.resources.MachineGroup(
 machine_type="c2d-highcpu-112",
 spot=True)
cloud_machine.start()
```

That’s it! That’s how simple it is to start a machine on the cloud with Inductiva - it only takes a 3-lines python script.
Run the script and you’ll have successfully started your first machine on the cloud!
You should now be able to see your active cloud machine in the <a href="https://console.inductiva.ai/machine-groups/active">Console</a>.
(If you’re struggling to create and run the python file, take a look at the <a href="quick-start-guide.html">Quick Start Guide</a>.)

## 3️⃣ Terminate your cloud machine

To wrap up, terminate your cloud machine, so that it doesn’t consume too many credits.

💰Reducing your costs is a priority for us, so Inductiva will automatically terminate any cloud machine after 3 minutes of idle time  (i.e. of not running any simulation task or executing associated services like downloading / uploading).
Nevertheless, you should always terminate the machine after the task finishes. You can do that either from the Console (click on your active machine in the list to enter its details, then click on the red button up top and confirm to shut it down); or by adding to the following instruction to your python script:

```python
cloud_machine.terminate()
```
