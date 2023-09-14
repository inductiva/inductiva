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

To test this scenario we have available one example of a protein PDB file - [download it here](inductiva/resources/alanine.pdb).

```python
 import inductiva

inductiva.api_key = "YOUR_API_KEY"

 # Initialize the scenario
 lysozyme = inductiva.molecules.examples.load_lysozyme()
 scenario = inductiva.molecules.ProteinSolvation(
     protein_pdb = lysozyme,
     temperature = 300)

 # Run a simulation
 task = scenario.simulate(simulation_time_ns = 0.5)

 # Get the simulation output on your local machine.
 output = task.get_output()
 ```
The last code line downloads the files necessary to render and post-process our simulation.

## Output and Post-processing 
The simulation output folder contains trajectory data spanning the simulation duration and the protein's topology. These files are necessary to generate an interactive trajectory visualization using [NGLview](https://github.com/nglviewer/nglview). Users have the flexibility to select the protein's representation type, including options to display or hide the backbone, and they can also specify the criteria for visualization selection. You can refer to the [documentation](https://nglviewer.org/ngl/api/manual/usage/selection-language.html) for guidance on creating custom selections according to your preferences.

```python
 output = task.get_output()
 output.render_interactive(representation="ribbon", add_backbone=True, selection="protein")
 ```
![protein_solvation](https://github.com/inductiva/inductiva/assets/114397668/e9dfe7d9-9c2b-4be4-90c3-0d68a2a7b9fd)

Users also have the option to plot the RMSF (Root Mean Square Fluctuation) over the simulation for each residue within the protein structure. [RMSF](https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html) measures the extent to which a structure deviates from its average configuration over time, offering insights into the mobility of specific protein residues. 

```python
 rmsf_values = output.plot_rmsf_per_residue()
 ```
<p align="center">
  <img src="https://github.com/inductiva/inductiva/assets/114397668/811bc191-6944-4b77-a0e2-0ba99b679246" alt="Centered Image" width="350" height="250">
</p>

Furthermore, you have the capability to visualize attributes for each amino acid, including metrics like RMSF (Root Mean Square Fluctuation), or any other properties that you compute independently.

```python
output.render_attribute_per_residue(rmsf_values)
 ```
<p align="center">
  <img src="https://github.com/inductiva/inductiva/assets/114397668/e9bbe012-030b-4e70-8774-fc2ff7bfb8e4" alt="Centered Image" width="350" height="250">
