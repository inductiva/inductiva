"""Class to run commands on an arbitrary image."""
from typing import Optional

from inductiva import types, tasks, simulators

class MpiConfig:
    """Class to configure MPI settings for the simulator."""

    def __init__(self, version: str = "",np: int = None, use_hwthread_cpus: bool = True):
        """Initialize the MPI configuration.
        Args:
            version: The version of the MPI library to use.
            np: Number of processes.
            use_hwthread_cpus: Use hyperthreading.
        """
        self.np = np
        self.version = version
        self.use_hwthread_cpus = use_hwthread_cpus

class Command:
    """Class to represent a command to run on an arbitrary image."""

    def __init__(self, command: str, mpi_config: Optional[MpiConfig] = None):
        """Initialize the Command class.
        Args:
            command: The command to run.
        """
        self.command = command
        self.mpi_config = mpi_config

@simulators.simulator.mpi_enabled
class ArbitraryImage(simulators.Simulator):
    """Class to run commands on an arbitrary image."""

    def __init__(self):
        """Initialize the ArbitraryImage class.
        Point to the API method to run a simulation.
        """
        super().__init__()
        self.api_method_name = "arbitrary.arbitrary_commands.run_simulation"

    def run(self,
            input_dir: str,
            commands: list[Command],
            container_image: Optional[str],
            storage_dir: Optional[str] = "",
            extra_metadata: Optional[dict] = None,
            on: Optional[types.ComputationalResources] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run.
            container_image: The container image to use for the simulation.
                Example: container_image="docker://inductiva/kutu:xbeach_v1.23"
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Parent directory for storing simulation
                               results.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           container_image=container_image,
                           **kwargs)

# import inductiva

# mg1 = inductiva.resources.MachineGroup(machine_type="c2d-standard-32")

# arb = inductiva.simulators.ArbitraryImage()

# config = MpiConfig(np=32, use_hwthread_cpus=True)

# command1 = Command("ls")
# command2 = Command("pwd", mpi_config=config)

# # Run on machine group with mpi
# #apptainer command 1
# #mpirun -np .... apptainer command 2
# arb.run(input_dir="folder1",
#         commmands=[command1, command2],
#         container_image="docker://inductiva/img1",
#         on= mg1)


# # Run on machine group without mpi
# #apptainer command 1
# #apptainer command 1
# arb.run(input_dir="folder1",
#         commmands=[command1, command1],
#         container_image="docker://inductiva/img1",
#         on= mg1)

# config = MpiConfig(np=32, use_hwthread_cpus=True, version="1.10.6")

# command3 = Command("ls", mpi_config=config)
# command4 = Command("pwd", mpi_config=config)

# # Run on MPICluster
# #setup the mpi version on all machines
# #mpirun -hosts file1 -np X .... apptainer command 3
# #mpirun -hosts file1 -np X .... apptainer command 4
# arb.run(input_dir="folder1",
#         commmands=[command3, command4],
#         container_image="docker://inductiva/img1",
#         on= mg1)