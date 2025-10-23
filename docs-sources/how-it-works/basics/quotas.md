# Quotas

Quotas help ensure fair access to powerful computing and storage for all users, so that users always get reliable performance when they launch their simulations.

You're free to use the Inductiva API as an individual user (individual plan) or join an
organization (enterprise or academia) - visit the [Pricing webpage](https://inductiva.ai/pricing) for more information.

| **Quota** | **Unit** | __Scope(*)__ | **Description** | **Individual** | **Enterprise** | **Academia** |
|-------|------|-------|-------------|----------|------------|------------|
| Maximum number of VCPUs | vcpu | global | Total number of VCPUs across all running machine instances plus the number of VCPUs of the instance to be requested must not exceed the quota limit | 720 | 10000 | 10000 |
| Maximum price per hour across all instances(**) | USD | global | Accumulated price per hour of all active machine instances, plus the price of the instance to be requested, must not exceed the quota | 10 | 1000 | 1000 |
| Maximum simultaneous instances | instance | global | Maximum number of machine instances running simultaneously at any moment | 10 | 100 | 100 |
| Maximum disk size | GB | instance | Maximum size of the disk that can be assigned to each individual machine in a machine group | 1000 | 5000 | 5000 |

**(*)NOTE:** _global_ quotas are applied to the user account and will encompass all
machines and tasks submitted by the user.
_Instance_ quotas are applied to each machine and task and, therefore,
are only applicable to each item individually.

__(**)NOTE:__  The Maximum price per hour for the active machine instances is not affected by the [Task Orchestration Fee](https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost#inductiva-s-task-orchestration-fee). I.e. this fee applied per simulation does not count for this quota's calculation.


## How to check quotas

The web Console's Dashboard shows the user's quotas, and their usage in real-time.
Also, there's a [CLI user command](https://inductiva.ai/guides/api-functions/cli/user)

If any of these quotas establish a limit for what you can achieve with Inductiva API, please reach out to us directly or ask for support via Discord ([join here](https://discord.com/invite/p9tjqBhuZ5) if you’re new, or [head to the server](https://discord.com/channels/1389190271723638804/1389235177456402502) if you’re already a member), so we can better understand your needs.

```{banner_small}
:origin: how_it_works_quotas
```
