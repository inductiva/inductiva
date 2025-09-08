# Setting Up the Simulation Environment

## Prerequisites
To follow this tutorial, either as a hands-on walkthrough or as a reference for your own implementation, please begin by downloading the simulation files [here]().

The files are organized into three folders to support both pre- and post-processing:

- `inputFiles_template` – contains the FEM model definition and analysis scripts
- `outputFiles` – the directory where the analysis results will be saved
- `Records` – contains the 30 acceleration time-histories required for the IDAs

A record duration file, `records_duration.txt`, is also included and will be read during the execution of the simulation script.

## Code Overview
The Python code provided below is designed to run OpenSees simulations on Inductiva’s cloud infrastructure. It is fully adaptable and can be customized for virtually any OpenSees use case.

In this example, we demonstrate how to run **50 Incremental Dynamic Analyses (IDAs) in parallel**, showcasing a typical high-throughput simulation workflow using the EESD OpenSees distribution.

The code is divided into six key sections:
1. Allocating the Cloud Machine Group
2. Setting Up the Project and Configuring the Simulator
3. Defining Analysis Parameters
4. Initializing the Simulation Loop and `TemplateManager`
5. Running the Simulations and Assigning Metadata
6. Monitoring Progress and Downloading Results

We’ll guide you through each of these sections to help you understand both the structure and flexibility of the workflow.

Copy the code and save it as a Python file - your simulation script - for execution after reviewing the other sections.

<code>