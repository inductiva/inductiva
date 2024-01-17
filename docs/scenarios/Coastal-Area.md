# Coastal Area Scenario

This scenario simulates waves propagating in a coastal area represented by a
bathymetric profile (i.e., the map of depths of the sea bottom as a function
of spatial coordinates x, y). Waves are injected into the domain from one of the
boundaries of the simulation with a given amplitude and period. The simulation
is performed using the [SWASH](https://swash.sourceforge.io/) simulator.

**Note:**  Some of the following tools require the installation of extra dependencies for the coastal package with
`pip install --upgrade "inductiva[coastal_extra]"`.

### Example 1: Simulating waves on a synthetic bathymetry

You can follow this example to simulate and visualize the wave propagation in a bathymetry
that represents a Coastal Area. To avoid having to download any data for now, we will start by simulating
on synthetic bathymetry, which is procedurally generated using the API.

```python
import inductiva

bathymetry = inductiva.coastal.Bathymetry.from_random_depths(x_range=(0,100),
          y_range=(0,100),
          x_num=100,
          y_num=100,
          max_depth=1,
          random_seed=12,
          initial_roughness=1,
          roughness_factor=0.5,
          percentile_above_water=20)

scenario = inductiva.coastal.CoastalArea(bathymetry=bathymetry,
                                             wave_source_location="W",
                                             wave_amplitude=0.1,
                                             wave_period=5.5)

task = scenario.simulate(simulation_time=80, output_time_step=1)

output = task.get_output()

output.render(movie_path="movie_path.mp4", fps=5)
```

Let's go in-depth to each of the steps in the example above.

We start by generating a random bathymetry using the [Diamond-Square algorithm](https://en.wikipedia.org/wiki/Diamond-square_algorithm).
The bathymetry has the following parameters that can be tuned by the user:
- `x_range`: The range of x values, in meters;
- `y_range`: The range of y values, in meters;
- `x_num`: Number of grid points in the x direction;
- `y_num`: Number of grid points in the y direction;
- `max_depth`: Maximum depth value, in meters;
- `initial_roughness`: Initial roughness value, in meters. Controls the
  initial range of randomness of the Diamond-Square algorithm;
- `roughness_factor`: Roughness factor. Must be between 0 and 1.
  Controls the rate at which the range of randomness of the
  Diamond-Square algorithm decreases over iterations. The higher the
  roughness_factor, the rougher the bathymetry;
- `percentile_above_water`: Percentile of the depths that must be above
  standard sea level, i.e., have a negative depth (or positive height) value.
  Must be between 0 and 100;
- `random_seed`: Since the Diamond-Square is a stochastic algorithm,
  we can specify a random seed to use.

The bathymetry object can be further plotted with:

```python
bathymetry.plot(show=True)
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/main/assets/bathymetry/bathymetry_random.png" alt="Raw bathymetry" width="550" height="450">
</p>

After the generation, the bathymetry is ready to be used to initialize the simulation scenario.
The user needs only to further define the following parameters:
- `water_level`: The water level, in meters. By default is 0 (standard sea-level);
- `wave_source_location`: The location of the wave source. Supported
  locations are N (north), S (south), E (east), and W (west),
  corresponding to the upper, lower, right and left boundaries of
  the simulation domain, respectively;
- `wave_amplitude`: The amplitude of the wave, in meters.
- `wave_period`: The period of the wave, in seconds.

Once the scenario is created, you can simulate by specifying the various
time parameters of the `scenario.simulate()` method:
- `simulation_time`: Total simulation time, in seconds;
- `time_step`: The system evolves in discrete time steps given by this argument;
- `output_time_step`: The system's information is written to the output files
(that contain the water level and velocities). This is not done for all time steps,
but only in `output_time_step` intervals.

To visualize the results, we can generate and save a movie of the simulation.

```python
output.render(movie_path = "movie_path.mp4", fps=5)
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/main/assets/media/random_coastal_area.gif" alt="Coastal area simulation" width="550" height="450">
</p>

### Example 2: Simulating waves on a real bathymetry

Now let's perform our simulations with real data!

Start by loading a bathymetry. This can be done e.g. by reading an
[ASCII XYZ file](https://emodnet.ec.europa.eu/sites/emodnet.ec.europa.eu/files/public/20171127_DTM_exchange_format_specification_v1.6.pdf),
a standard format for bathymetry data (for example, see e.g. the [European Marine Observation and Data Network](https://emodnet.ec.europa.eu/geoviewer/#!/)).
Here we shall use the bathymetry data from Algarve.

```python
import inductiva

bathymetry_url = "https://downloads.emodnet-bathymetry.eu/high_resolution/590_HR_Lidar_Algarve.emo.zip"

bathymetry_path = inductiva.utils.files.download_from_url(bathymetry_url)

bathymetry = inductiva.coastal.Bathymetry.from_ascii_xyz_file(
    ascii_xyz_file_path = bathymetry_path)
```

This bathymetry contains the depth data in 1683395 points, all across the entire Algarve
region, that spans over 120kms width. To visualize we can specify a resolution for the plot
(in this case 200m x 200m), which downsizes the bathymetry depths to a size more
manageable to the plot.

```python
bathymetry.plot(x_resolution=200, y_resolution=200, show=True)
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/main/inductiva/coastal/assets/media/bathymetry.png" alt="Algarve bath" width="550" height="450">
</p>

This raw bathymetry data spans an area too large and sparsely sampled to
be used in a simulation. In this case, the user can crop the bathymetry to a
smaller coastal area of interest:

```python
bathymetry = bathymetry.crop(x_range=(51000, 52000), y_range=(12150, 13000))
bathymetry.plot(show=True)
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/main/inductiva/coastal/assets/media/bathymetry_cropped.png" alt="Algarve bath" width="550" height="450">
</p>

To be used in a simulation, the bathymetry data must be interpolated to a
regular grid with custom resolution, in this case 5m x 5m. In case some
positions are missing in the bathymetry data (see e.g. the lower left corner in
the figure above), the user can choose to fill them with a constant depth or the
nearest depth value:

```python
bathymetry = bathymetry.to_uniform_grid(x_resolution=5,
                                        y_resolution=5,
                                        fill_value="nearest")
```

When interpolated to a regular grid, the bathymetry is ready to be used in a
simulation scenario. In this case, we set the wave source location to the south
boundary of the domain (i.e. the lower y boundary), with a wave amplitude of 5m
and a wave period of 5.5s:

```python
scenario = inductiva.coastal.CoastalArea(bathymetry=bathymetry,
                                         wave_source_location="S",
                                         wave_amplitude=5,
                                         wave_period=5.5)
```

In case the bathymetry is not defined on a regular grid when creating a
scenario, interpolation to a regular grid is automatically attempted with a
default resolution of 2m x 2m and no pre-configured strategy to fill depths at
positions missing in the bathymetry.

Once the scenario is created, you can simulate as follows:

```python
task = scenario.simulate(simulation_time=100, output_time_step=1)

output = task.get_output()
```

The user can specify the total simulation time and the time step between outputs
(all in seconds).

To visualize the results:

```python
output.render()
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/main/inductiva/coastal/assets/media/coastal_area.gif" alt="Algarve bath" width="550" height="450">
</p>

