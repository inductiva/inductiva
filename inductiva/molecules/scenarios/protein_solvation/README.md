# Protein Solvation Scenario

This scenario models the dynamics of a protein solvated in water inside a cubic box, powered by the GROMACs simulator. The system is modelled with the motion equations for each individual atoms according to the prescription of a molecular force field approximation.

To initialize the scenario, the user can define the following parameters:
- `protein_pdb` : PDB file containing the initial conformation of the protein;
- `temperature`: Temperature of the system in K.

Before simulation, files undergo automatic pre-processing in our machines to strip away all non-protein components of the system, such 
as ligands, non-standard residues, or water molecules, leaving only the protein.
We tested the scenario with structures predicted by the [AlphaFold2 model](https://alphafold.ebi.ac.uk) and structures from the [RCSB Protein Data Bank](https://www.rcsb.org), but you can provide any other well-formatted PDB file as input, as long as its contains at least one protein chain with all its atoms coordinates. 



The user can further define the following simulation parameters:
- `simulation_time_ns`: The simulation time, in nanoseconds. Default value is 10.
- `output_timestep_ps`: The time interval between recorded snapshots of the atom's positions in the trajectory output file, in picoseconds. Default value is 1.
- `integrator`: The integrator to use for the simulation.
    Options:
                - "md" ([Molecular Dynamics](https://manual.gromacs.org/documentation/2019/reference-manual/algorithms/molecular-dynamics.html)): Accurate leap-frog algorithm for
                integrating Newton's equations of motion.
                - "sd" ([Steepest Descent](https://manual.gromacs.org/current/reference-manual/algorithms/energy-minimization.html)): Stochastic dynamics integrator with
                leap-frog scheme.
                - "bd" ([Brownian Dynamics](https://manual.gromacs.org/documentation/2021.2/reference-manual/algorithms/brownian-dynamics.html)): Euler integrator for Brownian or
                position Langevin dynamics.
    Default value if "md".
- `n_steps_min`: Number of steps for energy minimization. Default value is 5000.

Moreover, the hardware and interaction are configured with the usual general parameters - `machine_group`, `run_async`, `n_cores`.
 Launching a simulation returns a task object, which can be used to verify the status of the simulation, get the simulation outputs and access post-processing tools. See more in [Tasks](inductiva/README.md).

### Example

```python
 import inductiva

 # Initialize the scenario
 scenario = inductiva.molecules.scenarios.ProteinSolvation(
     protein_pdb = "alanine.pdb"
     temperature = 300)

 # Run a simulation
 task = scenario.simulate(simulation_time_ns = 20)

 # Get the simulation output on your local machine.
 output = task.get_output()
 ```
