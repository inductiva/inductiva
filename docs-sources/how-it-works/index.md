# How It Works

At Inductiva, we want to provide you with the knowledge and tools you need to scale your simulations.
We've prepared a collection of helpful guides to support you navigating Inductiva's platform
with ease.

This section gives you a detailed overview of the inner workings of the Inductiva platform: from 
how simulations are structured to how your data is managed and scaled. It's designed to help you 
understand the building blocks so you can make the most out of Inductiva when running your simulations.

📘 What’s Inside:
* 🧱 [Building Blocks](building-blocks/index) – Explore the core structure of Inductiva—how Projects, 
Tasks, Machine Groups, and cloud Buckets connect. This foundational overview shows how simulations are 
organized and orchestrated, and how data flows within the platform.
* 🚀 [Get Started](get-started/index) – A practical guide to installing Inductiva and understand the 
onboarding flow. 
* 💰 [Costs and Quotas](basics/index) – Understand how credits translate into compute and storage usage, 
and see the quotas that define your resource limits. 
* 🧩 [Tasks](tasks/index) – Unpack everything related to simulation tasks from submission to status 
tracking and result outputs. Learn how tasks are managed and where to find key metadata to improve 
your workflow over time.
* ⚙️ [Computational Resources](machines/index) – A deep dive into how Machine Groups are defined, 
scaled, and benchmarked. Learn how to choose beteween the resources made available by our API, and 
understand when each option provides the best value for your specific use case.
* ☁️ [Cloud Storage](cloud-storage/index) – Get clarity on how Inductiva handles simulation data, 
where task outputs live, how to access and manage them, and why storage choices impact your credit 
usage. Learn best practices for organizing and cleaning up data.


💡 Why It’s Useful:
✓ Understand the platform’s architecture – Learn how Inductiva’s core components fit together so 
you can navigate and use them effectively.
✓ Launch and manage simulations with confidence – Get clear guidance from setup to output, whether 
you're running your first task or scaling up a complex project.
✓ Optimize performance and costs – Make informed decisions about machines, storage, and quotas to 
get faster results and better credit efficiency.
✓ Stay in control of your data and compute – Learn how simulations are executed, where your results 
are stored, and how to monitor or clean up your resources.
✓ Make the most of Inductiva’s flexibility – Whether you use built-in simulators or bring your own, 
you’ll understand how to tailor the platform to your workflow.


Explore the foundations and use them to build more powerful workflows: from setup to scaling, 
tuning to cleaning!


```{banner}
:origin: how_it_works_index
```

```{toctree}
---
caption: The Building Blocks
maxdepth: 2
hidden: true
---

Overview <building-blocks/index>
building-blocks/interfaces
building-blocks/configuring-simulators

```

```{toctree}
---
caption: Get Started  
maxdepth: 3
hidden: true
---

🚀 Install Inductiva API in 2 steps <get-started/install-guide>
⏩ Quick-Start guide <get-started/quick-start-guide>
📌 Pick a cloud machine for your simulation <get-started/pick-cloud-machine>
✈️ Start your first cloud machine with Inductiva <get-started/start-first-machine>
🛠️ Troubleshoot installation <get-started/troubleshooting>
🗑️ Uninstallation guide <get-started/uninstall_inductiva>

```

```{toctree}
---
caption: Costs and Quotas
maxdepth: 2
hidden: true
---

💲 How much compute time does 5$US provide <basics/compute-5usd>
💰 How much does a simulation cost in Inductiva? <basics/how-much-does-it-cost>
🔒 Inductiva Quotas <basics/quotas>
💥 How many cores can be used? <basics/how-many-cores>

```

```{toctree}
---
caption: Tasks
maxdepth: 3
hidden: true
---
Overview <tasks/index>
tasks/tasks
tasks/tasks-execution
tasks/tasks-lifecycle
tasks/manage_and_retrieve_results
```

```{toctree}
---
caption: Computational Resources
maxdepth: 3
hidden: true
---
Overview <machines/index>
machines/shared-dedicated-resources
machines/manage_computational_resources
machines/computational-infrastructure
machines/spot-machines
```

```{toctree}
---
caption: Cloud Storage
maxdepth: 2
hidden: true
---
Cloud Storage <cloud-storage/cloud-storage>
```
