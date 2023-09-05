# Protein Solvation Scenario

This scenario models the dynamics of a protein solvated in water inside a virtual cubic box. The system is modelled with the motion equations for each individual atoms according to the prescription of a molecular force field approximation.

To initilaize the scenario, the user can define the following parameters:
-- Protein .pdb file: structural .pdb file of the protein to simulate;
-- Temperature of the system in K.

Now, the user is ready to simulate the protein dynamics. Note that the scenario currently supports only protein .pdb files, like the ones available in AlphaFold. RCSB's protein files can be provided to this scenario, but all the non-proteic parts of the system like ligands, non-standard residues or water molecules will be stripped before simulation.

The user can further define the following simulation parameters:
-- simulator: The simulator to be used for the simulation. Currently, only `GROMACS` is supported.
-- simulation_time_ns: The simulation time, in nanoseconds. Default value is 10.
-- output_timestep_ps: The difference between the recorded snapshots of the atom's positions in the trajectory output file, in picoseconds. Default value is 1.
-- integrator: The integrator to use for the simulation.
    Options:
                - "md" (Molecular Dynamics): Accurate leap-frog algorithm for
                integrating Newton's equations of motion.
                - "sd" (Steepest Descent): Stochastic dynamics integrator with
                leap-frog scheme.
                - "bd" (Brownian Dynamics): Euler integrator for Brownian or
                position Langevin dynamics.
    Default value if "md".
-- n_steps_min: Number of steps for energy minimization. Default value is 5000.

Moreover, the hardware and interaction are configured with the usual general parameters - `machine_group`, `run_async`, `n_cores`.
 Launching a simulation returns a task object, which can be used to verify the status of the simulation, get the simulation outputs and access instantly post-processing tools. See more in [Tasks](inductiva/tasks/README.md).

### Example

```python
 import inductiva

 # Initialize the scenario
 scenario = inductiva.fluids.scenarios.ProteinSolvation(
     protein_pdb = "alanine.pdb"
     temperature = 300)

 # Run a simulation
 task = scenario.simulate(simulation_time_ns = 20)

 # Get the simulation output on your local machine.
 output = task.get_output()
 ```
