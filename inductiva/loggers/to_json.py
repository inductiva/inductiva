"""
BenchmarkLogger class responsible for writing the benchmark data to a json file.
"""
from typing import NamedTuple
import warnings
import json
import os
from datetime import datetime

import inductiva
from inductiva import types
from inductiva.loggers.benchmark_logger import BenchmarkLogger
from inductiva.client.model.task import Task
from inductiva.resources.machines_base import BaseMachineGroup



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
        }

        #create json file
        with open(self._path, "w", encoding="utf-8") as file:
            json.dump(data, file)

    def _update_json(self, data):
        """This function will update the json file with the new data.
        Replacing the old data with the new one.
        Args:
            data: Data to be updated in the json file."""
        with open(self._path, "r", encoding="utf-8") as file:
            old_data = json.load(file)

        #combine both dicst
        data = {**old_data, **data}

        with open(self._path, "w", encoding="utf-8") as file:
            json.dump(data, file)
    
    def _add_to_list(self, data, key):
        """This function will add data to a list in the json file.
        Args:
            data: Data to be added to the list.
            key: Key of the list in the json file."""
        with open(self._path, "r", encoding="utf-8") as file:
            old_data = json.load(file)

        #check if key exists
        if key not in old_data:
            old_data[key] = []

        old_data[key].append(data)

        with open(self._path, "w", encoding="utf-8") as file:
            json.dump(old_data, file)

    def log_benchmark(self, benchmark):
        """Logs benchmark.

        Args:
            benchmark: Name and session of the benchmark.
                benchmark = {
                    "name": str,
                    "session": int
                }
        """
        data = {
            "name": benchmark.name,
            "session": datetime.now().timestamp()
        }

        self._update_json(data)

    def log_resource(self, resource: BaseMachineGroup):
        """Logs resource to database.
        
        Args:
            resource: Resource to log.

        """
        data = {
        
            "machine_type": resource.name,
            "cost": resource.estimate_cloud_cost(),
            "machine_count": resource.num_machines,
            "provider": resource.provider,
            "id": resource.id
        }
        self._add_to_list(data, "resources")


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
        data = {
            "hash": test_case_instance.hash,
            "input_dir": test_case_instance.input_dir,
            "metadata": metadata
        }
        self._add_to_list(data, "testcases")
        pass

    def log_task(self, task: Task, test_case_instance: NamedTuple,
                 resource: BaseMachineGroup):
        """Logs information about a task.
        
        Args:
            task: Task to log.
            test_case: Test case associated with said task.
            resource: Resource that will be running the task.
        """
        data={
            "task_id": task.id,
            "testcase_hash": test_case_instance.hash,
            "resource_id": resource.id,
            "submit_time": datetime.now()
        }
        self._add_to_list(data, "tasks")

