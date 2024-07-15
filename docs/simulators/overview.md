# Overview

Inductiva API has available several ready-to-use open-source simulators. Users 
who are familiar with the simulators can easily start running simulations with 
their previously prepared simulation configuration files, and directly scale
them up on [last-generation hardware](https://tutorials.staging.inductiva.ai/intro_to_api/computational-infrastructure.html#available-computational-resources) 
avaialble on the cloud via Inductiva API.

The simulators available in the current version of the API (0.7) are:
- [AMR-Wind](../simulators/AmrWind.md)
- [CaNS](../simulators/CaNS.md)
- [DualSPHysics](../simulators/DualSPHysics.md)
- [FDS](../simulators/FDS.md)
- [GROMACS](../simulators/GROMACS.md)
- [OpenFOAM](../simulators/OpenFOAM.md)
- [OpenFAST](../simulators/OpenFAST.md)
- [Reef3D](../simulators/Reef3D.md)
- [SCHISM](../simulators/SCHISM.md)
- [SPlisHSPlasH](../simulators/SPlisHSPlasH.md)
- [SWAN](../simulators/SWAN.md)
- [SWASH](../simulators/SWASH.md)
- [XBeach](../simulators/XBeach.md)
- [NWChem](../simulators/NWChem.md)

## Typical Usage Patterns
As described in our [Tutorial](https://tutorials.inductiva.ai/intro_to_api/configuring-simulators.html)
typical usage patterns involve:
- [Single executable simulator;](https://tutorials.inductiva.ai/intro_to_api/configuring-simulators.html#the-simple-cases)
- [Pre-defined executables running in sequence;](https://tutorials.inductiva.ai/intro_to_api/configuring-simulators.html#a-slightly-more-complex-case)
- [Choose the executables as wished to run in sequence.](https://tutorials.inductiva.ai/intro_to_api/configuring-simulators.html#running-long-simulation-pipelines)

Check the documentation of each simulator to learn about the specifics of how
to configure them and launch your simulations via the API.

## Your favourite open-source simulator is not in the list above?
If you have any questions or suggestions about other open-source simulation
packages that you would like to see available via Inductiva API, please open an
issue on our [GitHub repository](https://github.com/inductiva/inductiva/issues),
or contact us via [support@inductiva.ai](mailto:support@inductiva.ai).
    
