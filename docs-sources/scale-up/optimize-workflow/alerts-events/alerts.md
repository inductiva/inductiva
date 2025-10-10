# 🔔 Alerts & Events

Inductiva provides a flexible alert system to keep users informed about important events in their simulations and account activity. Alerts can notify of issues that require attention, or confirm that expected events have occurred.
There are two types of alerts:
- Inductiva Alerts – predefined by Inductiva to warn users about events that could disrupt your work.
- Simulation Observer Events – user-defined alerts that notify users of specific events in their simulations.

Depending on the alert type, notifications may be sent by email, shown in the web Console, or both. The table below lists all Inductiva alerts, their purpose, and the communication channels through which they can be received.

## Inductiva Alerts

**Subject** | **Alert** | **Description** | **Communication channel** |
|---|---|---|---|
Tasks | Task success | Task completes successfully | Console notification |
Tasks | Task failed | Task fails or encounters errors | Console notification |
Tasks | Task output stalled | Task output stops updating | Console notification |
Tasks | Task preempted | Spot instances are preempted | Email + Console notification |
Machine Groups | Task runner error | Error in the component running a task | Console notification |
Machine Groups | Machine Group terminated: Low credits | Machine Group is terminated due to insufficient credits | Console notification |
Machine Groups | Quota exceeded | One of the [Inductiva Quotas](https://inductiva.ai/guides/how-it-works/basics/quotas) was exceeded | Console notification |
Machine Groups | VM preempted | Spot VM is reclaimed by the cloud provider | Email |
Machine Groups | Machine group pending start | Cloud provider has not allocated the requested resources. Pending request for over one hour | Email |
Credits | Credits exhausted | Credit balance reaches zero | Email + Console notification |
Credits | Credits below threshold | Current balance reaches 15% of the last top-up | Email + Console notification |
Credits | Credits below storage cost | Insufficient balance for monthly storage cost | Email + Console notification |
Credits | Top-up offer or referral | Credits added as a result of a campaign, referral, or other | Email |
System | Achievement unlocked | You unlocked new achievements | Console notification |
