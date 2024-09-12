# GROMACS

GROMACS is a versatile simulator to perform molecular dynamics simulations. It 
is primarily designed for biochemical molecules like proteins, lipids and
nucleic acids that have a lot of complicated bonded interactions, but since
GROMACS is extremely fast at calculating the nonbonded interactions (that
usually dominate simulations) many groups are also using it for research on
non-biological systems, e.g. polymers and fluid dynamics.

A single simulation of GROMACS via Inductiva API can comprise several steps - 
e.g., preparing the molecules, minimizing the energy of the system, running the
simulation and post-processing. Hence, to configure a simulation of GROMACS the
user may require several files. Moreover, GROMACS has specific commands to run
certain tasks that already use the files in your input directory. 

## Example

This example runs a simple case of the molecular dynamics of water molecules
inside a small box.

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-4",
    num_machines=1,
    data_disk_gb=10)
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "gromacs-input-example.zip", unzip=True)

commands = [
    "gmx solvate -cs tip4p -box 2.3 -o conf.gro -p topol.top",
    "gmx grompp -f energy_minimization.mdp -o min.tpr -pp min.top -po min.mdp -c conf.gro -p topol.top",
    "gmx mdrun -s min.tpr -o min.trr -c min.gro -e min.edr -g min.log",
    "gmx grompp -f positions_decorrelation.mdp -o decorr.tpr -pp decorr.top -po decorr.mdp -c min.gro",
    "gmx mdrun -s decorr.tpr -o decorr.trr -x  -c decorr.gro -e decorr.edr -g decorr.log",
    "gmx grompp -f simulation.mdp -o eql.tpr -pp eql.top -po eql.mdp -c decorr.gro",
    "gmx mdrun -s eql.tpr -o eql.trr -x trajectory.xtc -c eql.gro -e eql.edr -g eql.log",
]

gromacs = inductiva.simulators.GROMACS()

task = gromacs.run(input_dir=input_dir,
                   commands=commands,
                   on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
```
