**Find answers to commonly asked questions about COAWST.**

<br>

# FAQ

## How do I run a Delft3D simulation with FLOW-WAVE coupling?
To run a Delft3D simulation with FLOW-WAVE coupling using the Inductiva API,
you’ll need to run both the FLOW and WAVE components together. This involves
creating a shell script (e.g., `run_sim.sh`) that starts `d_hydro.exe` with
`mpirun` in the background and launches `wave.exe` in the foreground, allowing
them to interact during the simulation.

Here’s a basic example of such a script:

```bash
#!/bin/bash
flowexedir=$D3D_HOME/$ARCH/flow2d3d/bin
waveexedir=$D3D_HOME/$ARCH/wave/bin
swanexedir=$D3D_HOME/$ARCH/swan/bin
swanbatdir=$D3D_HOME/$ARCH/swan/scripts

# Run
# Start FLOW
mpirun -np 16 $flowexedir/d_hydro.exe config_d_hydro.xml &

# Start WAVE
$waveexedir/wave.exe r17.mdw 1
```

For a detailed walkthrough of the setup and execution process, refer to the full
this [tutorial](https://inductiva.ai/guides/delft3d/flow-wave-coupling).

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
