# Configuring OpenTelemac Simulations
OpenTelemac simulations using the Inductiva API are executed in a **single step**. Configuration is managed through a `.cas` steering file, along with additional files such as `.slf`, `.cli`, or `.bc`.

## The `.cas` File
The `.cas` file is the cornerstone of an OpenTelemac simulation, directing the solver's operations. For example, in our [Run Your First Simulation](https://inductiva.ai/guides/opentelemac/quick-start) tutorial, the `t2d_malpasset-fine.cas` file is used to configure the Malpasset Dam Break simulation.

### Key aspects controlled by the `.cas` file

- **Input/Output**: Defines geometry, boundary conditions, and result filenames.
- **Simulation Setup**: Specifies initial conditions, total simulation time (4000s), and a variable time step based on the Courant condition.
- **Physical Properties**: Includes bottom friction (Manning’s coefficient), turbulence modeling, and velocity diffusivity.
- **Numerical Scheme**: Implements the Saint-Venant equations with a first-order kinetic finite volume scheme, suitable for handling shocks and dry-bed wetting (enabled with `TIDAL FLATS = YES`).

### Example Lines from the `.cas` file

```
GEOMETRY FILE            = geo_malpasset-large.slf
BOUNDARY CONDITIONS FILE = geo_malpasset-large.cli
RESULTS FILE             = r2d_malpasset-fine.slf
EQUATIONS                = 'SAINT-VENANT FV'
FINITE VOLUME SCHEME     = 1
TIDAL FLATS              = YES
```


