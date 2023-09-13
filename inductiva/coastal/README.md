# Coastal dynamics

Inductiva API supports simulation of coastal dynamics via two simulators:
[SWASH](https://swash.sourceforge.io/) and [XBeach](https://oss.deltares.nl/web/xbeach/).
These simulators solve the [shallow-water equations](https://en.wikipedia.org/wiki/Shallow_water_equations),
a reduced version of the Navier-Stokes equations valid for small depth-to-length
ratios (i.e., shallow water).

You can invoke each of these simulators directly to execute your own simulation
scripts, or you can make use of the `CoastalArea` scenario, a high-level
interface for simulating wave conditions on bathymetries you can load.

## Invoking SWASH or XBeach directly

SWASH or XBeach can be invoked directly to run simulation input files you have
created. For example, to run a SWASH simulation which input files are stored in
the `simulation_input_dir` directory and whose main configuration file is
`INPUT.sws`, you can do:

```python
import inductiva

simulator = inductiva.simulators.SWASH()

task = simulator.run(
    input_dir="simulation_input_dir",
    sim_config_filename="INPUT.sws",
)

output = task.get_output()
```

Note that this runs the simulation in parallel, using 4 cores.

To run an XBeach simulation, the interface is the same, except that you need to
use the `XBeach` simulator class:

```python
simulator = inductiva.simulators.XBeach()
```

## Coastal area scenario

This scenario simulates waves propagating in a coastal area represented by a
bathymetric profile (i.e., the map of depths of the sea bottom as a function
of spatial coordinates x, y). Waves are injected into the domain from one of the
boundaries of the simulation with a given amplitude and period. The simulation
is performed using the [SWASH](https://swash.sourceforge.io/) simulator.

### Example

Start by creating or loading a bathymetry. This can be done e.g. by reading an
[ASCII XYZ file](https://emodnet.ec.europa.eu/sites/emodnet.ec.europa.eu/files/public/20171127_DTM_exchange_format_specification_v1.6.pdf),
a standard format for bathymetry data (for examples, see e.g. the [European Marine Observation and Data Network](https://emodnet.ec.europa.eu/geoviewer/#!/)).
Here we shall use the bathymetry data from Algarve available [here](https://sextant.ifremer.fr/record/SDN_CPRD_590_HR_Lidar_Algarve/)
as an example.

```python
import inductiva

bathymetry = inductiva.coastal.Bathymety.from_ascii_xyz_file(
    ascii_xyz_file_path="590_HR_Lidar_Algarve.emo")
```

The bathymetry can be inspected by plotting it with a given resolution, in this
case 200 m x 200 m:

```python
bathymetry.plot(x_resolution=200, y_resolution=200)
```

![Raw bathymetry.](resources/media/bathymetry.png)

Raw bathymetry data will often span an area too large and sparsely sampled to
be used in a simulation. In this case, the user can crop the bathymetry to a
smaller area of interest:

```python
bathymetry = bathymetry.crop(x_range=(51000, 52000), y_range=(12150, 13000))
bathymetry.plot()
```

![Cropped bathymetry.](resources/media/bathymetry_cropped.png)

To be used in a simulation, the bathymetry data must be interpolated to a
regular grid with custom resolution, in this case 5 m x 5 m. In case some
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
boundary of the domain (i.e. the lower y boundary), with a wave amplitude of 5 m
and a wave period of 5.5 s:

```python
scenario = inductiva.coastal.CoastalArea(bathymetry=bathymetry,
                                         wave_source_location="S",
                                         wave_amplitude=5,
                                         wave_period=5.5)
```

In case the bathymetry is not defined on a regular grid when creating a
scenario, interpolation to a regular grid is automatically attempted with a
default resolution of 2 m x 2 m and no pre-configured strategy to fill depths at
positions missing in the bathymetry.

Once the scenario is created, run the simulation as follows:

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

![Coastal area simulation.](resources/media/coastal_area.gif)