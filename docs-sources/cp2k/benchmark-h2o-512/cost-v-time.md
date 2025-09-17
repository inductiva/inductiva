# Cost vs Time

This plot compares the **execution cost** against the **simulation runtime**, helping identify the **sweet spot** between **speed** and **cost efficiency**.

```{raw} html
:file: ../_static/cost_time_graph.html
```

> **Note:**
> The plot shows a significant increase in both cost and simulation time when moving from the `c4-highcpu-144` to the `c4-highcpu-192`. This may be due to the fact that C4 machines can be equipped with one of two CPU types—[Granite Rapid or Emerald Rapid](https://cloud.google.com/compute/docs/general-purpose-machines#c4_series). The slower simulation is likely running on a different CPU architecture, which could explain the large performance difference.

```{banner_small}
:origin: cp2k_benchmark-512_cost-v-time
```