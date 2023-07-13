
[![Python package](https://github.com/inductiva/inductiva/actions/workflows/python-package.yml/badge.svg)](https://github.com/inductiva/inductiva/actions/workflows/python-package.yml)

![linkedin_header](https://user-images.githubusercontent.com/104431973/231184851-0ce34289-593e-4832-aaa2-9aae652113f5.jpg)

## Large scale simulations made simple

Inductiva API provides dozens of open-source physical simulation packages from your laptop. Users can start simulating right away, with no hardware setup issues and no software configuration headaches. Additionally, we provide a transparent way to scale your simulations to the next-level with one-line of code.

One can **run** simulations in two ways:
- **High-level Simulation:** We provide pre-built scenarios that define a physical model (e.g., the dam break scenario in fluid dynamics). Users are allowed to define some parameters of the scenarios, choose the method of simulation with a few parameters and accelerate their simulations with the best hardware available for them. 

Example of how to run the dam break scenario:
```python
from inductiva import fluids

scenario = fluids.scenarios.DamBreak(dimensions=(1., 0.3, 0.3))

output = scenario.simulate()
video = output.render()
```

- **Low-level Simulation:** Users familiar with the simulators can submit their previously prepared simulation configuration files. In this way, users can immediatelly take advantage of performant hardware to speed up their simulations, without having to change any of their existing simulation scripts. Moreover, users may want to execute more than just one simulation and via API they can do it in a couple of lines of code. 

Example of how to run a low-level simulation:
```python
from inductiva import fluids

simulator = fluids.simulators.DualSPHysics()

output_dir = simulator.run(input_dir="FlowCylinder",
                           sim_config_filename="CaseFlowCylinder_Re200_Def.xml",
                           output_dir="Flow",
                           device="gpu")
```

Find more examples of simulations at the [tutorials section](https://github.com/inductiva/inductiva/tree/main/demos).

Our goal is to provide researchers and engineers with an easy and fast way to scale their simulations and explore various designs. 


## Simulators

The simulators we provide are all open-source and have their own dedicated documentation.

Currently, we have available the following simulators:
- [SPlisHSPlasH](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH)
- [DualSPHysics](https://github.com/DualSPHysics/DualSPHysics)
- [OpenFOAM](https://www.openfoam.com/)
- [SWASH](https://swash.sourceforge.io/)
- [xBeach](https://oss.deltares.nl/web/xbeach/)
- [GROMACS](https://www.gromacs.org/)

If you would like other simulators to be added, contact us at [simulations@inductiva.ai](mailto:simulations@inductiva.ai).

## Installation

It is super simple to start using the API if you are already familiar with Python package management.
One just needs to do
```
pip install inductiva
```

and your are good to go! You are ready to start [exploring our tutorial notebooks](https://github.com/inductiva/inductiva/tree/main/demos).

Notice that, in a local computer you just need to do this once. When opening the tutorials in Google Colab, you may need to re-install this
between notebooks.

## API access tokens

Please request your demo API token via the following Google Form (we will reply to you by email):

[Request API token](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform)

