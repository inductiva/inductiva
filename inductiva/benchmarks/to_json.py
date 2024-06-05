"""
BenchmarkLogger class responsible for writing the benchmark data to a json file.
"""
from typing import Dict, NamedTuple, Union
import warnings
import json
import os

import inductiva
from inductiva import types
from inductiva.benchmarks.benchmark_logger import BenchmarkLogger
from inductiva.client.model.task import Task
from inductiva.resources.machines_base import BaseMachineGroup
from inductiva.utils import files


class ToJson(BenchmarkLogger):
    """Class responsible for writing the benchmark data to a json file."""

    def __init__(self, path: types.Path = None):
        self._path = path

    def create_json(self, json_name: str):
        output_dir = inductiva.get_output_dir()

        if os.path.isdir(output_dir):
            warnings.warn("Path already exists, files may be overwritten.")
        else:
            os.makedirs(output_dir)

        self._path = f"{str(output_dir)}/{json_name}.json"

        data = {
            "benchmark_name": json_name,
        }

        #create json file
        with open(self._path, 'w', encoding='utf-8') as file:
            json.dump(data, file)

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

    def log_resource(self, resource: BaseMachineGroup):
        """Logs resource to database.
        
        Args:
            resource: Resource to log.

        """
        pass

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

    def log_task(self, task: Task, test_case_instance: NamedTuple,
                 resource: BaseMachineGroup):
        """Logs information about a task.
        
        Args:
            task: Task to log.
            test_case: Test case associated with said task.
            resource: Resource that will be running the task.
        """
        pass
