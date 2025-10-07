# Part 1: Case Setup

## 1. Defining the Damping Function
We begin by defining a **damping function** that calculates the coefficients α (mass-proportional damping) and β (stiffness-proportional damping) according to the Rayleigh method, which are essential for modeling energy dissipation in structural systems.

This function takes the periods of two vibration modes and a target damping ratio as input. It returns coefficients that ensure consistent damping behavior across a specified frequency range. These coefficients are fundamental in OpenSees simulations for realistic modeling of energy dissipation under dynamic loading.

```python
from os.path import join
import numpy as np
import inductiva
import math

# Defining the Damping Function
def damping(Ti, Tj, ksi):
    """
    Calculate damping coefficients alpha and beta.

    Parameters:
        Ti (float): Period of the first mode.
        Tj (float): Period of the second mode.
        ksi (float): Damping ratio.

    Returns:
        tuple: (alpha, beta) damping coefficients.
    """
    fi = 1 / Ti
    fj = 1 / Tj

    wi = 2 * math.pi * fi
    wj = 2 * math.pi * fj

    alpha = ksi * 2 * wi * wj / (wi + wj)
    beta = ksi * 2 / (wi + wj)

    return alpha, beta
```

## 2. Project and File Path Configuration
To manage and organize simulation data, we group everything into a **Project**. We also define the input file paths required for the analyses.

Don't forget to update the `tutorial_folder` variable to point to the path of the downloaded IDA-at-scale folder.

```python
# Organize all data into a Project
project_name = "B2_3P_Tutorial"

# ---- Path's ----
# Path to the Tutorial Folder
tutorial_folder = r"/Path/to/IDA-at-scale"

# Name of the main tcl file
tcl_file = "Prototype_b2_3p_batch.tcl"

# Path to the template folder inside the tutorial
input_files_template = join(tutorial_folder, "inputFiles_template")
# Path to the simulation files after render inside the tutorial
simulation_files = join(tutorial_folder, "simulation_files")
# Path to the records file inside the tutorial
records_duration_file = join(tutorial_folder, "records_duration.txt")
# Load records into memory
records_duration = np.loadtxt(records_duration_file, delimiter=' ')
# ---- Paths ----
```

> Learn more about **Projects** [here](https://inductiva.ai/guides/scale-up/projects/index).

## 3. Defining Analysis Parameters
We now set the parameters for the nonlinear dynamic analyses. The IDA will run 30 earthquake ground motion records, each scaled to 10 intensity levels using predefined EQfactor values.

We also define numerical damping, required for non-linear dynamic analyses. Rayleigh damping is specified in OpenSees via the Rayleigh command, which needs two parameters, alpha and beta, based on two fundamental frequencies of the structure. In this example, the Rayleigh damping is calibrated based on the 1<sup>st</sup> and 6<sup>th</sup> frequencies of the building. A damping ratio of 5% is assumed.

```python
# ---- Analysis variables ----
# 30 bidirectional acceleration time-series to perform IDA
analysis_range = range(1,31)

# Scaling factor for the earthquake time-series
EQfactor_values = [
    0.05, 0.10, 0.15,
    0.20, 0.25, 0.30,
    0.35, 0.40, 0.45,
    0.50]

# Computing the parameters for a 5% fraction
# damping based on the 1st and 6th frequencies
damping_percentage = 0.05
Ti = 0.09970895596567399 # 1st vibration period
Tj = 0.04970502055069268 # 6th vibration period

alpha, beta = damping(Ti, Tj, damping_percentage)
# ---- Analysis variables ----
```

With all case-related settings complete, let's explore how to run the 300 simulations on the cloud using Inductiva.