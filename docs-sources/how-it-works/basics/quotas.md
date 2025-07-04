## Quotas

You're free to use the Inductiva API as an individual user (starter plan) or join an
organization (business or academia) - visit the Pricing webpage for more information .

| Quota | unit | scope* | Description | Starter | Business | Academia |
|-------|------|-------|-------------|----------|------------|------------|
| Maximum number of VCPUs | vcpu | global | Total number of VCPUs across all running machine instances plus the number of VCPUs of the instance to be requested must not exceed the quota limit | 300 | 10000 | 10000 |
| Maximum price per hour across all instances | USD | global | Accumulated price per hour of all active machine instances, plus the price of the instance to be requested, must not exceed the quota | 10 | 1000 | 1000 |
| Maximum simultaneous instances | instance | global | Maximum number of machine instances running simultaneously at any moment | 10 | 100 | 100 |
| Maximum disk size | GB | instance | Maximum size of the disk that can be assigned to each individual machine in a machine group | 1000 | 5000 | 5000 |

***NOTE:** _global_ quotas are applied to the user account and will encompass all
machines and tasks submitted by the user.
_Instance_ quotas are applied to each machine and task and, therefore,
are only applicable to each item individually.

## How to monitor your account details

Inductiva's CLI provides an easy way to monitor your account details and get
information about your tier, credits and current quotas. Simply use
the command `inductiva user info`


If any of these quotas establish a limit for what you can achieve with Inductiva
API, please [reach out to us](mailto:support@inductiva.ai) and we can better
understand your needs.
