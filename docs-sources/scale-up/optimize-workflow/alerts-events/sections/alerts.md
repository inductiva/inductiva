# Inductiva Alerts

Inductiva Alerts are predefined by Inductiva to warn users about events that could disrupt their work.


**Subject** | **Alert** | **Description** | **Communication channel** |
|---|---|---|---|
Tasks | Task success | Task completes successfully | Console notification |
Tasks | Task failed | Task fails or encounters errors | Console notification |
Tasks | Task output stalled | Task output stops updating | Console notification |
Tasks | Task preempted | Spot instances are preempted | Email + Console notification |
Machine Groups | Task runner error | Error in the component running a task | Console notification |
Machine Groups | Machine Group terminated: Low credits | Machine Group is terminated due to insufficient credits | Console notification |
Machine Groups | Quota exceeded | One of the [Inductiva Quotas](https://inductiva.ai/guides/how-it-works/basics/quotas){:target="_blank"} was exceeded | Console notification |
Machine Groups | VM preempted | Spot VM is reclaimed by the cloud provider | Email |
Machine Groups | Machine group pending start | Cloud provider has not allocated the requested resources. Pending request for over one hour | Email |
Credits | Credits exhausted | Credit balance reaches zero | Email + Console notification |
Credits | Credits below threshold | Current balance reaches 15% of the last top-up | Email + Console notification |
Credits | Credits below storage cost | Insufficient balance for monthly storage cost | Email + Console notification |
Credits | Top-up offer or referral | Credits added as a result of a campaign, referral, or other | Email |
System | Achievement unlocked | You unlocked new achievements | Console notification |


Are we missing an important event that you'd like to be alerted of? Contact us on [Inductiva’s Discord](https://discord.com/invite/p9tjqBhuZ5)


```{banner_small}
:origin: recipes_delete_failed
```