# Benchmark

## Three-Dimensional Currents
This benchmark focuses on the _S1 simulation_ as outlined in the paper "[Modeled Three-Dimensional Currents and Eddies on an Alongshore-Variable Barred Beach](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020JC016899)", authored by _Christine M. Baker, Melissa Moulton, Britt Raubenheimer, Steve Elgar, and Nirnimesh Kumar_ (2021). To better suit the purposes of this benchmark, we limited the simulation duration to 10 and 100 seconds to ensure that each simulation completes within approximately 5 and 28 minutes, respectively, when running on the slowest machines.


```{note}
_The original input files can be found [here](https://zenodo.org/records/4091612)._
```



### Software Versions
The software versions used in this benchmark were the following:

| Component              | Version                               |
|------------------------|---------------------------------------|
| SWASH                  | v9.01A                                |
| Clang                  | 16.0.0 (clang-1600.0.26.6)            |
| kernel                 | 24.5.0                                |


### Execution times comparison
Faster execution times are better. Units in seconds.

```{raw} html
:file: time_graph.html
```

### Cost vs Time
The following plot shows how the execution times and costs vary with the number of cores for each machine family.

```{raw} html
:file: cost_time_graph.html
```