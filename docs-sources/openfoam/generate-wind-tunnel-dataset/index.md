# Generate a Wind Tunnel Simulation Dataset
In this tutorial, we will demonstrate how to use the Inductiva API to efficiently generate an OpenFOAM dataset 
in parallel by varying the inlet wind speed. This approach allows rapid creation of multiple simulation cases 
with different wind conditions, facilitating the analysis of flow behavior and performance sensitivity over 
a controlled range of wind speeds. By leveraging cloud resources, users can run these simulations in parallel 
and accelerate the generation of CFD data for further study or integration into downstream workflows.

> 📺 **Prefer video?**  
> This guide is also available as a [webinar replay](../webinars/openfoam-cfd-dataset.md) where we walk through running **OpenFOAM on Inductiva** step by step.  
> [_Watch it to see the process in action!_](../webinars/openfoam-cfd-dataset.md)

<p align="center"><img src="../_static/bike_pressure_field.png" alt="OpenFOAM motorBike visualization" width="700"></p>

To demonstrate this process, we will use the incompressible, steady-state simpleFoam solver along with 
the **motorBike example**, which is available in the OpenFOAM repository.

In this tutorial, we'll walk through how to:
- [Review the prerequisites for running the motorBike case](https://inductiva.ai/guides/openfoam/generate-wind-tunnel-dataset/sections/section1)
- [Run the simulation](https://inductiva.ai/guides/openfoam/generate-wind-tunnel-dataset/sections/section2)
- [Generalize the use case](https://inductiva.ai/guides/openfoam/generate-wind-tunnel-dataset/sections/section3)
- [Generate synthetic wind tunnel data at scale](https://inductiva.ai/guides/openfoam/generate-wind-tunnel-dataset/sections/section4)
- [Postprocessing with Inductiva](https://inductiva.ai/guides/openfoam/generate-wind-tunnel-dataset/sections/section5)
- [Results and Key Takeaways](https://inductiva.ai/guides/openfoam/generate-wind-tunnel-dataset/sections/section6)

```{banner_small}
:origin: openfoam
```


```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
sections/section5.md
sections/section6.md
```