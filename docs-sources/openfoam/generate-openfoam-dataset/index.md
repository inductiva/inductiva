# Generate a OpenFOAM dataset
[OpenFOAM](https://www.openfoam.com/) is a powerful open-source CFD toolbox widely used across engineering and scientific disciplines for simulating complex fluid dynamics problems. With capabilities extending from incompressible and compressible flows to heat transfer, turbulence modeling, multiphase flows, and more, OpenFOAM supports high-fidelity simulations in both research and industrial applications.

In this tutorial, we will demonstrate how to use the Inductiva API to efficiently generate an OpenFOAM dataset in parallel, varying the inlet wind speed. This approach enables rapid creation of multiple simulation cases with different wind conditions, making it easier to analyze flow behavior and performance sensitivity across a controlled range of wind speeds. By leveraging cloud resources, users can run these simulations in parallel and accelerate the generation of CFD data for further study or integration into downstream workflows.

<p align="center"><img src="../_static/bike_pressure_field.png" alt="OpenFOAM simulation visualization" width="700"></p>

To demonstrate this, we will use the incompressible, steady-state simpleFoam solver with the [motorBike example](https://develop.openfoam.com/Development/openfoam/-/tree/master/tutorials/incompressible/simpleFoam/motorBike), available on the [OpenFOAM  repository](https://develop.openfoam.com/Development/openfoam).

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
sections/section5.md
```
