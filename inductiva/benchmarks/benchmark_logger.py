"""Abstract class for benchmark loggers."""
from abc import ABC, abstractmethod
from collections import namedtuple
from typing import Dict, NamedTuple, Union

from inductiva.client.model.task import Task
from inductiva.resources.machines_base import BaseMachineGroup


class BenchmarkLogger(ABC):
    """Abstract class for benchmark loggers."""

    test_case_instance = namedtuple("TestCaseInstance",
                                    ["input_dir", "hash", "testcase", "id"])

    @abstractmethod
    def log_benchmark(self, benchmark: Dict[str, Union[int, str]]):
        """Logs benchmark.

        Args:
            benchmark: Name and session of the benchmark.
                benchmark = {
                    "name": str,
                    "session": int
                }
        """
        pass

    @abstractmethod
    def log_resource(self, resource: BaseMachineGroup):
        """Logs resource to database.
        
        Args:
            resource: Resource to log.

        """
        pass

    @abstractmethod
    def log_testcase(self, test_case_instance: NamedTuple, metadata: dict):
        """Logs a testcase.

        This method saves the task input id (unique value representing the
            task input files), the task name and the metadata used on said task.
        Args:
            test_case_instance: NamedTuple with the input_dir, hash, testcase
                and id. This hash is the unique identifier
                generated from the input_dir.
            metadata: Metadata to run the task.
        """
        pass

    @abstractmethod
    def log_task(self, task: Task, test_case_instance: NamedTuple,
                 resource: BaseMachineGroup):
        """Logs information about a task.
        
        Args:
            task: Task to log.
            test_case: Test case associated with said task.
            resource: Resource that will be running the task.
        """
        pass
