# How much does it cost to simulate in Inductiva?

Inductiva acts as an intermediary between the user and the cloud provider, so the cost of the resources are directly linked to their cost at the cloud provider (e.g. Google Cloud Provider) and is broken down into:

| Cost  | What is it? | When is it applied? | How's it calculated? | How much is it? (Approximately) |
| --------------- | ------------- | ------------- | ------------- | ------------- |
| Compute |  Cost of the cloud machine where the simulation is run | While the cloud machine is active, regardless if it's running a simulation or not | Price per hour (*) times Machine's Uptime | A core-hour costs approx. 0.01\$US in Spot mode (**); and between 0.04\$US and 0.07\$US in regular mode |
| Storage | Cost of keeping files in the user's cloud bucket | While there are files stored in the cloud bucket | Price per GB per month x How long it's been stored | Around 0.02\$US per GB per month |
| Data Export (Downloading) |  Cost of extracting data from the cloud provider, such as downloading it to a local computer or sending it to another cloud | Whenever the user exports data; OR or when shared data is exported by someone else. | Price per GB x Volume of data exported | Around 0.12\$US per GB |
| Data Import (Uploading) |  Cost of uploading data to the cloud provider | N/A | N/A | No cost |

(*) The price/hour of the machines depend on their attributes, namely the number of vCPUs and the amount of RAM per vCPU. The price/hour fora machine with the same attributes may still vary depending on the globe region where it's hosted. Inductiva grants access to hundreds of machines, and empowers the user to select the best option for your simulation use-case - see how to Pick your cloud machine.   

(**) **Spot** instances are a type of virtual machine offering at significantly reduced costs compared to regular instances, with the downside that they can be terminated at any time. This characteristic makes spot ideal for budget-conscious users when cost savings outweigh the risk of unexpected termination, namely if the task is expected to finish in a short timeframe.

## Inductiva's fee 
Inductiva's fee comes solely from the percentage applied when purchasing credits.
1. Users purchase credits whenever they want via Stripe;
2. We deduct a Platform fee (%) from each top-up, so your available credits will be the top-up amount minus the fee;
3. The Compute, Storage and Data Transfer costs are subtracted from the credits only while they are being used.

Check out the <a href="https://inductiva.ai/pricing">Pricing page</a> for more information about our price model.


## How much can you accomplish with 5\$US
How much compute time will 5\$ get you at Inductiva? Well, 5\$ can actually get you quite a lot compute time!
 
💡Consider for example **Inductiva’s top pick machine**, the c2d-highcpu-112. 
This 112-core machine with 224GM of RAM can cost you as little as 0.68\$US per hour in spot mode. So, 5\$US will let you use this machine for more than 6 hours, which is equivalent to about 700 core-hours. 


## How do you control your costs
**Inductiva's Cost transparency policy** guarantees full visibility on how your credits are being spent.   
Here's the level of granularity with which the costs are presented throughout Inductiva's platform:
+ Available machines - Price/hour of the machines available on Inductiva;
+ Started machines - Estimated cost based on the machine's price/hour and its uptime;
+ Terminated machines - Cost of the machine based on its price/hour and uptime;
+ Tasks - Estimated cost based on its duration and the cost of the machine where it was run;
+ Storage - Cost of keeping the files in the user's cloud bucket;
+ Monthly Accumulated Costs - split into Compute, Storage and Data Transfer.

