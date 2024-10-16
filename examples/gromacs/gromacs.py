"""GROMACS example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "gromacs-input-example.zip",
    unzip=True)

commands = [
    "gmx solvate -cs tip4p -box 2.3 -o conf.gro -p topol.top",
    ("gmx grompp -f energy_minimization.mdp -o min.tpr -pp min.top -po min.mdp "
     "-c conf.gro -p topol.top"),
    "gmx mdrun -s min.tpr -o min.trr -c min.gro -e min.edr -g min.log",
    ("gmx grompp -f positions_decorrelation.mdp -o decorr.tpr -pp decorr.top "
     "-po decorr.mdp -c min.gro"),
    ("gmx mdrun -s decorr.tpr -o decorr.trr -x  -c decorr.gro -e decorr.edr "
     "-g decorr.log"),
    ("gmx grompp -f simulation.mdp -o eql.tpr -pp eql.top -po eql.mdp "
     "-c decorr.gro"),
    ("gmx mdrun -s eql.tpr -o eql.trr -x trajectory.xtc -c eql.gro -e eql.edr "
     "-g eql.log"),
]

gromacs = inductiva.simulators.GROMACS()

task = gromacs.run(input_dir=input_dir, commands=commands, on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
