"""gprMax example.

This example demonstrates three ways to run GprMax:
1. Without MPI (single process)
2. With GprMax's built-in MPI (using -mpi flag)
3. With mpirun wrapper (using --mpi-no-spawn flag)
"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
gprmax = inductiva.simulators.GprMax(version="3.1.7")


# OPTION 1: Run without MPI (single process)
# -------------------------------------------------
commands_no_mpi = [
    "python -m gprMax user_models/cylinder_Bscan_2D.in -n 60",
]

task = gprmax.run(
    input_dir="/Path/to/gprmax-input-example",
    commands=commands_no_mpi,
    on=cloud_machine)


# OPTION 2: Run with GprMax's built-in MPI
# -------------------------------------------------
# Use this mode when you want GprMax to manage MPI internally.
# The -mpi flag tells GprMax how many processes to use.
commands_gprmax_mpi = [
    "python -m gprMax user_models/cylinder_Bscan_2D.in -n 60 -mpi 61",
]

task = gprmax.run(
    input_dir="/Path/to/gprmax-input-example",
    commands=commands_gprmax_mpi,
    use_gprmax_mpi=True,  # Enable GprMax's built-in MPI mode
    on=cloud_machine)


# OPTION 3: Run with mpirun wrapper
# -------------------------------------------------
# Use this mode when you want to use standard mpirun.
# The --mpi-no-spawn flag is required to prevent GprMax from spawning
# its own MPI processes.
commands_mpirun = [
    "python -m gprMax antenna_like_MALA_1200_fs.in -n 3 --mpi-no-spawn",
]

task = gprmax.run(
    input_dir="/Path/to/gprmax-input-example",
    commands=commands_mpirun,
    n_vcpus=4,  # Specify number of MPI processes
    use_gprmax_mpi=False,  # Use mpirun wrapper (default)
    on=cloud_machine)


# OPTION 4: Mixed workflow (simulation + post-processing)
# -------------------------------------------------
# Common pattern: Run simulation with MPI, then post-process without MPI.
# The MPI configuration is ONLY applied to GprMax simulation commands.
# Utility commands (tools.*) automatically run without MPI.
commands_mixed = [
    # Step 1: Run simulation with MPI (4 processes)
    "python -m gprMax user_models/cylinder_Bscan_2D.in -n 60 --mpi-no-spawn",
    # Step 2: Merge output files (NO MPI - runs on single process)
    "python -m tools.outputfiles_merge user_models/cylinder_Bscan_2D",
    # Step 3: Plot results (NO MPI - runs on single process)
    "python -m tools.plot_Bscan user_models/cylinder_Bscan_2D_merged.out Ez",
]

task = gprmax.run(
    input_dir="/Path/to/gprmax-input-example",
    commands=commands_mixed,
    n_vcpus=4,  # MPI only for the simulation command
    on=cloud_machine)


# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
