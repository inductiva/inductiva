# IDAs on Inductiva: Run 300 Simulations in Parallel

*This tutorial was written by* [Daniel Caicedo](mailto:caicedod93@gmail.com) *in collaboration with the* **Inductiva Team**

<br>

**Incremental Dynamic Analysis (IDA)** is a computational method widely used in Performance-
Based Earthquake Engineering (PBEE) to evaluate the seismic performance of structures. This 
method generates response curves across a range of earthquake intensity levels, enabling a
probabilistic assessment of structural behaviour and associated seismic risk.

To date, the applicability of this technique for the case of unreinforced masonry constructions has been limited due to the **complexity of numerical models** and the **large computational burden** associated.

The first challenge is addressed by developing efficient structural models in OpenSees based on three-
dimensional macroelements to capture the in-plane (IP) and out-of-plane (OOP) mechanisms,
while also incorporating the effects of non-linear floor-to-wall connections and wall-to-wall
interlocking. 

The second, and more significant, challenge is tackled by leveraging the Inductiva
platform, which enables the execution of large-scale simulations through high-performance
computing (HPC) resources.

This tutorial showcases an example involving a two-storey unreinforced masonry building,
representative of the pre-code masonry building stock in metropolitan area of Lisbon. All the data
required for the development of numerical models was gathered in the scope of the [STAND4HERITAGE](https://stand4heritage.org) project (New STANDards for seismic assessment of built cultural HERITAGE).

<image2 and caption - ask Caicedo>

## Prerequisites
Download the required files [here](<link to input files>) and place them in a folder called `SimulationFiles`. 

-> Are these really prerequisites?
- Output files folder to save the results of the analysis.
- Records folder, where the 30 acceleration time-histories required for the IDAs are
located.
- Record duration file to be read during the execution of the script.

## Run 50 Simulations in Parallel
This section demonstrates how to run 50 Incremental Dynamic Analyses (IDAs) in parallel using an elastic group of cloud machines. The process is divided into six main code sections.

### Section 1: Allocating the Cloud Machine Group
We begin by allocating an elastic machine group to run the simulations. This group is configured with 
a minimum and maximum number of active machines. Machines are dynamically turned on or off
as demanded, ensuring efficient resource utilization.

Since some simulations may take longer than others, this elastic setup prevents idle waiting. When 
a machine finishes its assigned tasks and no further simulations are queued, it automatically shuts 
down. This dynamic scaling optimizes cloud usage and helps reduce computational costs.

```python
```

-> I think we should explain the machine we're using 

### Section 2: Project Setup and Simulator Configuration
First, we specify the directory for the designated input files folder.

Next, we specify the OpenSees version to use. In this case, version 3.0.2, which includes the 3D macroelement, OrthotropicMembraneSection, and NoTensionSection3d commands, and supports the use of non-symmetric tangent stiffness in ND-zeroLength elements.

To ensure all simulation data is easily managed and to simplify post-processing and result downloads, it's best to associate all data with a project using the Inductiva Project class.

> Learn more about the Inductiva Project class [here](https://inductiva.ai/guides/scale-up/projects/projects).

```python
```

### Code Section 3: Computing Parameters for the Analysis
Next, we define the analysis range for the 30 records used in the IDA, along with the corresponding EQfactor values, which denote the level of intensity to scale the ground motion records.

We also define numerical damping, required for non-linear dynamic analyses. Rayleigh damping is specified in OpenSees via the Rayleigh command, which needs two parameters, alpha and beta, based on two fundamental frequencies of the structure. In this example, the Rayleigh damping is calibrated based on the 1<sup>st</sup> and 6<sup>th</sup> frequencies of the building. A damping ratio of 5% is assumed, typical for this type of analysis.

```python
```

### Code Section 4: Initializing the Loop and `TemplateManager`
This section starts the loops over `analysis_range` and `EQfactor_values`. Each record is scaled to the specified intensity level, meaning the input files are rewritten at each iteration with updated parameters.

To automate this process, the `TemplateManager` tool is used. It overwrites the input files (with a .jinja extension) by replacing placeholders ({{ }}) with the appropriate parameter values at each step.

> Learn more about the `TemplateManager` [here](https://website-staging.inductiva.ai/guides/scale-up/parallel-simulations/templating).

```python
```

### Code Section 5: Running Simulations and Defining Metadata
Within the loops, each simulation is executed by calling the batch file. Metadata is also defined for each run, capturing information essential for future post-processing.

```python
```

### Code Section 6: Monitoring Simulations and Downloading Results
In the final stage, we wait for all simulations to complete. Once finished, the output data is automatically downloaded for further analysis, and the cloud machine group is terminated.

```python
```

## Post-Processing and Analysis of Results
Once the simulation data is available, post-processing can be performed using tools of your choice—such as custom Python scripts.

In this tutorial, results are analyzed in terms of:
1. Capacity and fragility functions
2. Convergence of statistical moments (specifically, the median and dispersion)

Capacity curves are presented based on the 50th, 16th, and 84th percentiles of the IDA response, using the maximum global drift ratio as the engineering demand parameter (EDP). The 50th fractile represents the median capacity of the structure.

Next, the damage states of Immediate Occupancy (IO), Collapse Prevention (CP), and Near Collapse (NC) are identified. The IO limit is defined approximately at the end of the elastic branch. CP is identified when the slope of the capacity curve decreases to 20% of the initial elastic slope. Finally, NC is indicated by the appearance of the flatline.

<image2 and caption - ask Caicedo>

Fragility functions denote the likelihood of exceedance of a specified EDP value corresponding to each limit state, as shown in the figure below.

<image3 and caption - ask Caicedo>

The figures below illustrate the convergence of the median and dispersion of the fragility functions, based on the 30 observations for each limit state.

<image4/5 and caption - ask Caicedo>

## -> I think it's lacking a conclusion or an ending section