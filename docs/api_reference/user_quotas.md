# User Tiers, Capabilities and Quotas

Users of the Inductiva API are subject to a set of rules that define the
resources they can use and the limits they are subject to. These rules are
defined by the user's tier, the set of capabilities they have access to, and the
quotas that limit the amount of resources they can use.
In this document, we will explain how these rules are defined and how they
affect the user's experience with the Inductiva API.

## Tiers

The Inductiva API has four tiers: Freemium, Standard, Power-user, and Enterprise.
The following table shows the main differences between these tiers:

| Tier | Description |
|------|-------------|
| Freemium | The Freemium tier is the entry-level tier for the Inductiva API. It is free to use and provides access to a limited set of capabilities and resources. Members of this tier can only submit tasks to the shared queue, _i.e_, they cannot allocate and use a dedicated machine group for their tasks. |
| Standard | The Standard tier provides access to a wider set of capabilities and resources than the Freemium tier. Members can allocate dedicated machine groups but with limited capabilities and low quota limits. |
| Power-user | The Power-user tier extends on the capabilities and resources available to the Standard tier but with higher quota limits and access to a wider set of capabilities. |
| Enterprise | The Enterprise tier provides access to all capabilities and resources available in the Inductiva API. Members of this tier have the highest quota limits and can use all the features available in the API. |

## Capabilities


| Capability| Description | Freemium | Standard | Power-user | Enterprise |
|-----------|-------------|----------|----------|------------|------------|
|Allow tasks on "default machine group" | Ability to submit tasks to the shared queue - machine group shared among all users that runs tasks in order of submission. | ✅︎ | ✅︎ |✅︎ |✅︎ |
| Allow tasks on dedicated machine group | Ability to launch and allocate a machine group, to be dedicated to the the user’s tasks. | ❌ | ✅︎ | ✅︎ | ✅︎ |
| Allow the override of time to live for tasks | Ability to extend the maximum time that a task can stay running, delaying its automatic termination. | ❌ | ✅︎ | ✅︎ | ✅︎ |
| Allow running "non Kutu" containers | Ability to specify the use of third-party docker containers that are defined/built outside the context of the KUTU repository. | ❌ | ❌ | ❌ | ✅ |
| Allow use of the simulators - any version | Ability to run any version of all simulators available in Inductiva API | ✅︎ | ✅︎ | ✅︎ | ✅︎ |
| Allow the use of MPI cluster | Ability to allocate multiple machines configured as a single cluster for running parallelized simulations across multiple machines using MPI | ❌ | ❌ | ❌ | ✅︎ |
| Allow the use of elastic clusters | Ability to dynamically allocate a machine pool with flexible sizing | ❌ | ✅︎ | ✅︎ | ✅︎ |
| Allow the use of "spot" resources | Ability to request spot resources (typically machine groups) to run tasks. Spot resources are resources that can be preempted at any moment by the cloud provider at any moment | ❌ | ✅︎ | ✅︎ | ✅︎ |
| Allow the use of "on demand" resources | Ability to request dedicated resources (typically machine groups) to run tasks. | ❌ | ✅︎ | ✅︎ | ✅︎ |
| Automatically restart tasks interrupted by preemption | A task that was running in a “spot” instance that was preempted (taken back by cloud provider) will be  automatically resubmitted | ❌ | ❌ | ✅︎ | ✅︎ |
| xxxx | xxxx | ❌ | ✅︎ | ✅︎ | ✅︎ |

## Quotas

| Quota | unit | scope | Description | Freemium | Standard | Power-user | Enterprise |
|-------|------|-------|-------------|----------|----------|------------|------------|
| Maximum tasks per week | task | global |Total number of tasks ran in the last 7-days window, including the task to be submitted, must not exceed the quota limit| 30 | 300 | inf | inf |
| Maximum number of VCPUs | vcpu | global | Total number of VCPUs across all running machine instances plus the number of VCPUs of the instance to be requested must not exceed the quota limit | N/A | 160 | 1000 | inf |
| Maximum price per hour across all instances | USD | global | Accumulated price per hour of all active machine instances, plus the price of the instance to be requested, must not exceed the quota | 2 | 4 | 270 | inf |
| Maximum simultaneous instances | instance | global | Maximum number of machine instances running simultaneously at any moment | N/A | 40 | 100 | inf |
| Maximum time a machine group can stay idle before termination | minute | global | Maximum time that an allocated machine group is allowed to stay active without running a task | N/A | 60 | 120 | 240 |
| Maximum time a machine group can stay up before automatic termination | hour | instance | Maximum time that an allocated machine group is allowed to stay active (even if it’s running a task), after which it will be automatically terminated | N/A | 36 | 48 | inf |
| Maximum time a task can stay running before automatic termination | hour | instance | Maximum time that a task can stay running, after which it will be automatically terminated | 4 | 8 | 16 | inf |
| Maximum disk size | GB | instance | Maximum size of the disk that can be assigned to each individual machine in a machine group | N/A | 100 | 200 | inf |
| Maximum amount of RAM per VCPU | GB | instance | Maximum amount of RAM per individual VCPU that can be used. Even though RAM is not specifiable per se, this quota constrains the machine types that can be requested | N/A | 4 | 6 | 8 |


## FAQs:

**When are tiers, credits and quotas verified?**

> Every time the user tries to register a MG or submits a task.

***When are credits subtracted from the user’s current amount?***

> Every time a task reaches a final state, the cost of the task is calculated and stored in the database. If the task is running on a shared resource the cost must be substracted to the user balance. Every time a machine is stopped, the cost of the machine is calculated and stored in the database, the machine cost must be substracted to the user balance.

***How are tiers, capabilities and quotas linked?***

> Tier basically defines which of the API functionalities are available to users, and which limits (quotas) are the users subject to.
> There are multiple capabilities available in the Inductiva API, as well as several limits (quotas) on multiple computational aspects. 
> These capabilities and quotas may change when you enroll in a certain campaign, such as Genesis.
> Users can be enrolled in multiple campaigns, but they are assigned to just one Tier.

***What happens if a user tries accessing capabilities or resources that are not within their tier’s limits?***

> The task is not executed and an error message is printed to the CLI. 
> If it was caused by quota overpassing, the message explicitly informs which quota would be over the limit if the task would run, so that the user is able to act upon this information and reattempt to run the task.



## How to monitor your account details


Inductiva's CLI provides an easy way to monitor your account details, and get
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
 Maximum time a task can stay running in the default queue before automatic termination   16 hour
 Maximum disk size                                                                        2000 GB
 Maximum amount of RAM per VCPU                                                           6 GB

```

This information is also available programatically through the python client:

```python
import inductiva

inductiva.users.get_info() # <-- to get information about the tier and credits
inductiva.users.get_quotas() # <-- to get quotas

```






If any of these quotas establish a limit for what you can achieve with Inductiva
API, please [reach out to us](mailto:support@inductiva.ai) and we can better
understand your needs.
