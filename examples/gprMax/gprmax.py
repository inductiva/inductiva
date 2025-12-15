"""gprMax example.

GprMax supports multiple ways to run with MPI.

See: https://docs.gprmax.com/en/latest/openmp_mpi.html
"""
import inductiva
from inductiva.commands import Command, MPIConfig

# Allocate cloud machine
machine = inductiva.resources.MachineGroup(provider="GCP",
                                           machine_type="c3d-highcpu-180")

gprmax = inductiva.simulators.GprMax(version="3.1.7")

# ============================================================================
# EXAMPLE 1: No MPI (single process)
# ============================================================================
commands = ["python -m gprMax cylinder_Bscan_2D.in -n 60"]

task = gprmax.run(input_dir="./gprmax-files", commands=commands, on=machine)

# ============================================================================
# EXAMPLE 2: GprMax's built-in MPI
# ============================================================================
# GprMax spawns its own MPI processes using the -mpi flag
commands = ["python -m gprMax cylinder_Bscan_2D.in -n 60 -mpi 61"]

task = gprmax.run(input_dir="./gprmax-files", commands=commands, on=machine)

# ============================================================================
# EXAMPLE 3: Using mpirun wrapper
# ============================================================================
# Use --mpi-no-spawn to prevent GprMax from spawning its own processes
# Wrap the command with MPIConfig to control mpirun
commands = [
    Command("python -m gprMax antenna_like_MALA_1200_fs.in -n 3 --mpi-no-spawn",
            mpi_config=MPIConfig(version="4.1.6", np=4))
]

task = gprmax.run(input_dir="./gprmax-files", commands=commands, on=machine)

# ============================================================================
# EXAMPLE 4: Full workflow with post-processing
# ============================================================================
# Mix simulation and post-processing commands as needed
commands = [
    # Run simulation with GprMax's built-in MPI
    "python -m gprMax cylinder_Bscan_2D.in -n 30 -mpi 31",
    "python -m tools.outputfiles_merge cylinder_Bscan_2D",
]

task = gprmax.run(input_dir="./gprmax-files", commands=commands, on=machine)

task.wait()
machine.terminate()
