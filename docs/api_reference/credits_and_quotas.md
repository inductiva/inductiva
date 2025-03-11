# Credits and Quotas

All users receive 5 US$ in credits upon registration. As a user you may also
purchase credits via Stripe whenever you want.

Credits are consumed based on the user's resource usage. This includes factors
like computation time and the specific machine groups utilized. For example,
running a simulation on a high-memory machine group will consume more credits
than on a standard machine group.

Real-time tracking of credit usage allows users to track and optimize their
consumption efficiently and avoid unexpected shortfalls. Below we'll
describe the new tools that guarantee that users have full visibility over their
available credits. If a user has insufficient credits to launch a new resource,
a clear message will be displayed.

By understanding and managing their credits, users can optimize their resource
usage ensuring they get the most value from their allocated credits while
staying within their limits.

## Quotas

You're free to use the Inductiva API as an individual user or join an
organization (business or academic).

| Quota | unit | scope* | Description | Individual | Business | Academic |
|-------|------|-------|-------------|----------|------------|------------|
| Maximum number of VCPUs | vcpu | global | Total number of VCPUs across all running machine instances plus the number of VCPUs of the instance to be requested must not exceed the quota limit | 300 | 10000 | 10000 |
| Maximum price per hour across all instances | USD | global | Accumulated price per hour of all active machine instances, plus the price of the instance to be requested, must not exceed the quota | 10 | 1000 | 1000 |
| Maximum simultaneous instances | instance | global | Maximum number of machine instances running simultaneously at any moment | 10 | 100 | 100 |
| Maximum disk size | GB | instance | Maximum size of the disk that can be assigned to each individual machine in a machine group | 1000 | 5000 | 5000 |

***NOTE:** _global_ quotas are applied to the user account and will encompass all
machine groups and tasks submitted by the user.
_Instance_ quotas are applied to each machine group and task and, therefore,
are only applicable to each item individually.

## FAQs

**When are credits and quotas verified?**

> Every time the user tries to register a machine group or submits a task.

## How to monitor your account details

Inductiva's CLI provides an easy way to monitor your account details and get
information about your tier, credits and current quotas. Simply use
the command `inductiva user info`:

```bash
$ inductiva user info
Name: <name of the user here>
Email: <user e-mail here>
Username: <username here>

■ Credits: 1000.00 US$

■ Global User quotas
                                                                 CURRENT USAGE     MAX ALLOWED
 Maximum number of VCPUs                                         0 vcpu            1000 vcpu
 Maximum price per hour across all instances                     0 USD             270 USD
 Maximum simultaneous instances                                  0 instance        100 instance

■ Instance User quotas
                                                                                          MAX ALLOWED
 Maximum disk size                                                                        1000 GB

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
