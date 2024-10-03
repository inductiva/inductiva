# Tiers, Capabilities and Quotas

Each user is assigned to a specific tier, which comes with an associated number
of credits. The tier defines the functionalities accessible to the user. Each tier
is associated with quotas that establish limits on the computational resources a
user can utilize.

Credits are consumed based on the user's resource usage. This includes factors
like computation time and the specific machine groups utilized. For example,
running a simulation on a high-memory machine group will consume more credits
than on a standard machine group. Users spend their credits on the resources and
functionalities available within their tier. The credits continue to be consumed
until they are either exhausted or reach their expiration date (if applicable).

Real-time tracking of credit usage allows users to track and optimize their
consumption efficiently and avoid unexpected shortfalls. Below we'll
describe the new tools that guarantee that users have full visibility over their
current tier and available credits. If a user attempts to access resources or
functionalities that they are not allowed to due to their tier or have insufficient
credits for, a clear message will be displayed. This notification will inform
the user about the restriction and suggest possible actions, such as upgrading
their tier or adjusting their resource usage.

By understanding and managing their credits, users can optimize their resource usage
ensuring they get the most value from their allocated credits while staying within
their limits.

In this document, we will explain how these rules are defined and how they
affect the user's experience with the Inductiva API.

## Tiers

The Inductiva API has three tiers: **Standard**, **Power-user**, and **Enterprise**.
The following table shows the main differences between these tiers:

| Tier | Description |
|------|-------------|
| Standard | The Standard tier serves as the entry-level option for all members, allowing the allocation of dedicated machine groups. However, it comes with limited functionality and lower quota restrictions. |
| Power-user | The Power-user tier extends on the capabilities and resources available to the Standard tier but with higher quota limits and access to a wider set of capabilities. |
| Enterprise | The Enterprise tier provides access to all capabilities and resources available in the Inductiva API. Members of this tier have the highest quota limits and can use all the features available in the API. |

## Capabilities

Each tier comes with a set of capabilities that define the functionalities
available to the user. The following table shows the main capabilities available
in each tier:

| Capability| Description | Standard | Power-user | Enterprise |
|-----------|-------------|----------|------------|------------|
| Allow tasks on dedicated machine group | Ability to launch and allocate a machine group, to be dedicated to the the user’s tasks. | ✅︎ | ✅︎ | ✅︎ |
| Allow the override of time to live for tasks | Ability to extend the maximum time that a task can stay running, delaying its automatic termination. | ✅︎ | ✅︎ | ✅︎ |
| Allow running "non Kutu" containers | Ability to specify the use of third-party docker containers that are defined/built outside the context of the KUTU repository. | ❌ | ❌ | ✅ |
| Allow use of the simulators - any version | Ability to run any version of all [simulators available in Inductiva API](../simulators/overview.md) | ✅︎ | ✅︎ | ✅︎ |
| Allow the use of standard machine groups | Ability to allocate a machine pool with a fixed number of machines | ✅︎ | ✅︎ | ✅︎ |
| Allow the use of elastic clusters | Ability to dynamically allocate a machine pool with flexible sizing | ❌ | ✅︎ | ✅︎ |
| Allow the use of MPI cluster | Ability to allocate multiple machines configured as a single cluster for running parallelized simulations across multiple machines using MPI | ❌ | ❌ | ✅︎ |
| Allow the use of "spot" resources | Ability to request spot resources (typically machine groups) to run tasks. Spot resources are resources that can be preempted at any moment by the cloud provider at any moment | ✅︎ | ✅︎ | ✅︎ |
| Allow the use of "on demand" resources | Ability to request dedicated resources (typically machine groups) to run tasks. | ✅︎ | ✅︎ | ✅︎ |
| Automatically restart tasks interrupted by preemption | A task that was running in a “spot” instance that was preempted (taken back by cloud provider) will be automatically resubmitted | ❌ | ✅︎ | ✅︎ |

## Quotas

| Quota | unit | scope* | Description | Standard | Power-user | Enterprise |
|-------|------|-------|-------------|----------|------------|------------|
| Maximum tasks per week | task | global |Total number of tasks ran in the last 7-days window, including the task to be submitted, must not exceed the quota limit| 300 | inf | inf |
| Maximum number of VCPUs | vcpu | global | Total number of VCPUs across all running machine instances plus the number of VCPUs of the instance to be requested must not exceed the quota limit | 160 | 1000 | inf |
| Maximum price per hour across all instances | USD | global | Accumulated price per hour of all active machine instances, plus the price of the instance to be requested, must not exceed the quota | 4 | 270 | inf |
| Maximum simultaneous instances | instance | global | Maximum number of machine instances running simultaneously at any moment | 40 | 100 | inf |
| Maximum time a machine group can stay idle before termination | minute | global | Maximum time that an allocated machine group is allowed to stay active without running a task | 60 | 120 | 240 |
| Maximum time a machine group can stay up before automatic termination | hour | instance | Maximum time that an allocated machine group is allowed to stay active (even if it’s running a task), after which it will be automatically terminated | 36 | 48 | inf |
| Maximum time a task can stay running before automatic termination | hour | instance | Maximum time that a task can stay running, after which it will be automatically terminated | 8 | 16 | inf |
| Maximum disk size | GB | instance | Maximum size of the disk that can be assigned to each individual machine in a machine group | 100 | 200 | inf |
| Maximum amount of RAM per VCPU | GB | instance | Maximum amount of RAM per individual VCPU that can be used. Even though RAM is not specifiable per se, this quota constrains the machine types that can be requested | 4 | 6 | 8 |

***NOTE:** _global_ quotas are applied to the user account and will encompass all
machine groups and tasks submitted by the user.
_Instance_ quotas are applied to each machine group and task and, therefore,
are only applicable to each item individually.

## FAQs

**When are tiers, credits and quotas verified?**

> Every time the user tries to register a machine group or submits a task.

***When are credits subtracted from the user’s current amount?***

> Every time a task reaches a final state, the cost of the task is calculated and stored in the database. If the task is running on a shared resource the cost must be substracted to the user balance. Every time a machine is stopped, the cost of the machine is calculated and stored in the database, the machine cost must be substracted to the user balance.

***How are tiers, capabilities and quotas linked?***

> Tier basically defines which of the API functionalities are available to users, and which limits (quotas) are the users subject to.
> There are multiple capabilities available in the Inductiva API, as well as several limits (quotas) on multiple computational aspects. 
> These capabilities and quotas may change when you enroll in a certain campaign, such as Genesis.
> Users can be enrolled in multiple campaigns, but they are assigned to just one Tier.

***What happens if a user tries accessing capabilities or resources that are not within their tier’s limits?***

> The task is not executed and an error message is printed to the CLI.
> If it was caused by quota overpassing, the message explicitly informs which quota would be over the limit if the task would run, so that the user is able to act upon this information and reattempt to run the task after changing the parameters.

## How to monitor your account details

Inductiva's CLI provides an easy way to monitor your account details and get
information about your tier, credits and current quotas. Simply use
the command `inductiva user info`:

```bash
$ inductiva user info
Name: <name of the user here>
Email: <user e-mail here>
Username: <username here>

■ Tier: Power-User

■ Credits

  Power-User (tier)               0.00
  pioneer (campaign)          10000.00
  ------------------------------------
  Total                       10000.00

■ Campaigns

 NAME      ENROLLMENT DATE     EXPIRY DATE        AVAILABLE CREDITS     INITIAL CREDITS
 pioneer   2024-02-06 11:40    2024-07-31 01:00   10000                 10000

■ Global User quotas
                                                                 CURRENT USAGE     MAX ALLOWED
 Maximum tasks per week                                          0 task            N/A
 Maximum number of VCPUs                                         0 vcpu            1000 vcpu
 Maximum price per hour across all instances                     0 USD             270 USD
 Maximum simultaneous instances                                  0 instance        100 instance
 Maximum time a machine group can stay idle before termination   N/A               120 minute

■ Instance User quotas
                                                                                          MAX ALLOWED
 Maximum time a machine group can stay up before automatic termination                    48 hour
 Maximum disk size                                                                        2000 GB
 Maximum amount of RAM per VCPU                                                           6 GB

```

This information is also available programmatically through the Python client:

```python
import inductiva

inductiva.users.get_info() # <-- to get information about the tier and credits
inductiva.users.get_quotas() # <-- to get quotas

```


If any of these quotas establish a limit for what you can achieve with Inductiva
API, please [reach out to us](mailto:support@inductiva.ai) and we can better
understand your needs.

