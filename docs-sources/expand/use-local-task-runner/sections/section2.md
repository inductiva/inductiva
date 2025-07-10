# Step 2: Run a Simulation Locally
With the task-runner active, you can now run a simulation on your local machine. Below is an example using the GROMACS simulator:

```python

import inductiva

# Specify the machine group
machine_group = inductiva.resources.machine_groups.get_by_name(
    machine_name="<machine-group-name>")

# Download the input files
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "gromacs-input-example.zip",
    unzip=True
)

# Define the simulation commands
commands = [
    "gmx solvate -cs tip4p -box 2.3 -o conf.gro -p topol.top",
    "gmx grompp -f energy_minimization.mdp -o min.tpr -pp min.top -po min.mdp -c conf.gro -p topol.top",
    "gmx mdrun -s min.tpr -o min.trr -c min.gro -e min.edr -g min.log",
    "gmx grompp -f positions_decorrelation.mdp -o decorr.tpr -pp decorr.top -po decorr.mdp -c min.gro",
    "gmx mdrun -s decorr.tpr -o decorr.trr -x -c decorr.gro -e decorr.edr -g decorr.log",
    "gmx grompp -f simulation.mdp -o eql.tpr -pp eql.top -po eql.mdp -c decorr.gro",
    "gmx mdrun -s eql.tpr -o eql.trr -x trajectory.xtc -c eql.gro -e eql.edr -g eql.log",
]

# Initialize the GROMACS simulator and run the simulation
gromacs = inductiva.simulators.GROMACS()
task = gromacs.run(
    input_dir=input_dir,
    commands=commands,
    on=machine_group)

# Wait for the task to complete
task.wait()
```

This script shows how to configure and run a molecular dynamics simulation using GROMACS locally through the Inductiva platform.