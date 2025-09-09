## Results and Key Takeaways
Once the simulation data is available, post-processing can be performed using tools of your choice, such as custom Python scripts.

In this tutorial, results are analyzed in terms of:
1. Capacity and fragility functions
2. Convergence of statistical moments (specifically, the median and dispersion)

Capacity curves are presented based on the 50th, 16th, and 84th percentiles of the IDA response, using the maximum global drift ratio as the engineering demand parameter (EDP). The 50th fractile represents the median capacity of the structure.

Next, the damage states of Immediate Occupancy (IO), Collapse Prevention (CP), and Near Collapse (NC) are identified. The IO limit is defined approximately at the end of the elastic branch. CP is identified when the slope of the capacity curve decreases to 20% of the initial elastic slope. Finally, NC is indicated by the appearance of the flatline.

<p align="center"><img src="../../_static/graph1.png" alt="Capacity curves with limit state observations" width="700"></p>

Fragility functions denote the likelihood of exceedance of a specified EDP value corresponding to each limit state, as shown in the figure below.

<p align="center"><img src="../../_static/graph2.png" alt="Fragility functions" width="700"></p>

The figures below illustrate the convergence of the median and dispersion of the fragility functions, based on the 30 observations for each limit state.

<p align="center"><img src="../../_static/graph3.png" alt="Convergence of fragility median" width="700"></p>

<p align="center"><img src="../../_static/graph4.png" alt="Convergence of fragility dispersion" width="700"></p>

Together, these results give a clear and thorough picture of the structure’s performance and reliability during seismic events.

Performing IDAs to this level of detail typically demands significant computational time. However, with Inductiva’s cloud platform, the process is completely transformed. What normally would take **over 150 hours** to run sequentially can be done in just **5 hours and 28 minutes** by running 200 `c2d-highcpu-2` machines in parallel.

This shows how Inductiva can dramatically speed up large-scale structural simulations, helping you get faster insights while using computational resources more efficiently.