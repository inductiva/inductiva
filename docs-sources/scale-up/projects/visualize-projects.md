# Visualize Projects

The Inductiva [Web Console](https://console.inductiva.ai/dashboard) provides a visual interface for managing projects, complementing the Python API and CLI commands.

## Project Creation

To create a new project in the console:
1. Click the _"Create Project"_ button in the [projects section](https://console.inductiva.ai/projects).

!["Create new Project in the Inductiva Console"](./_static/create_project.png)

2. In the dialog box that appears, enter a name for your project.
3. Click _"Create"_ to finalize. Your new project is now ready.

!["Save new Project in the Inductiva Console"](./_static/create_project_2.png)

## Project Visualization
Selecting a project takes you to its dedicated dashboard that gives you a real-time, at-a-glance overview of your simulations progress and resource usage.

### Project Dashboard Overview
The dashboard is split into two main areas: summary cards at the top for quick insights and detailed tables below.

#### Summary Cards

These cards highlight the most important metrics of your project:

- **Success Rate**: The percentage of tasks that completed successfully.
- **Average Task Duration**: The average runtime for a single task.
- **Active / Total Tasks**: The number of currently running tasks versus the total number of tasks submitted.
- **Active / Total Machine Groups**: The number of machine groups currently active versus the total number launched for the project.
- **Estimated Compute Cost**: Estimated computational cost in US dollars ($) for all tasks in the project.

#### Tables

For more granular information, the dashboard includes two tables:
- **Tasks Table**: Detailed breakdown of every task. You can check each task's `ID`, status (e.g., `Success`, `Failed`), duration, the simulator used, and its individual estimated cost.
- **Machine Groups Table**: Lists all machine groups launched for the project. It shows the machine type used (e.g., `c4-highcpu-4`), the total duration it was active, and its status (e.g., `Terminated`, `Running`).

!["Project Dashboard in the Inductiva Console"](./_static/project_dashboard.png)

### Accessing Task Details

To dive deeper into a specific simulation, click on its row in the **Tasks table**. This will take you to the Task Detail page:

## Project Management
### Renaming Project

### Moving Tasks

### Deleting Project