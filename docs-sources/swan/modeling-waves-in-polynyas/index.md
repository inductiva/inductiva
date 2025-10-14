# Modeling Waves in Antarctic Polynyas
In this tutorial, we demonstrate how to use the Inductiva API to run a series of SWAN simulations that replicate scenarios from the study: [Herman, A., & Bradtke, K. (2023). SWAN wave model simulations of the Terra Nova Bay Polynya](https://zenodo.org/records/8308164).

This example is particularly compelling because it walks through a **complete scientific modeling workflow** — from preparing real-world data, to running simulations, to validating the results. Simulating wave growth in Antarctic polynyas is a complex challenge due to the combination of strong offshore winds, rapidly forming frazil and grease ice, and sharp spatial gradients between open water and ice-covered regions. These conditions require precise input preparation and careful parameter tuning, as even small changes can significantly affect the simulated wave spectra and their agreement with satellite observations.

Using Inductiva, such studies can be greatly streamlined. Once your input data is prepared, you can run multiple simulations in parallel on high-performance cloud machines — with minimal scripting effort. This significantly reduces the time needed to complete large-scale simulation studies, helping you move faster from model to results.

In this tutorial, we’ll focus on two key simulation scenarios:
- **S0** — Open-water baseline (no ice effects). This is used to generate reference spectra for comparison with ice-affected simulations.
- **S2_f5** — Full ice scenario, where frazil and grease ice are included, and an ice-related dissipation term is activated with a specific parameter setting (power-law exponent = 5).

Let's get started!

```{banner_small}
:origin: swan
```

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
```
