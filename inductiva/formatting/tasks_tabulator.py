"""
Tabulator for lists of tasks.
"""
from .base_tabulator import BaseTabulator, Col
from ..utils import format_utils
from inductiva.client import models

from inductiva import tasks

emph_formatter = format_utils.get_ansi_formatter()


def status_formatter(status, unused_task):
    if status == models.TaskStatusCode.SUCCESS:
        return emph_formatter(status, format_utils.Emphasis.GREEN)
    elif status in tasks.Task.FAILED_STATUSES:
        return emph_formatter(status, format_utils.Emphasis.RED)
    return status


def datetime_formatter(dt, unused_task):
    return format_utils.datetime_formatter(dt)


def header_formatter(header: str):
    return emph_formatter(header.upper(), format_utils.Emphasis.BOLD)


def computetime_fmtr(_, task):
    execution_time = task.get_computation_time(cached=True)

    if execution_time is not None:
        execution_time = format_utils.seconds_formatter(execution_time)
        if task.info.computation_end_time is None:
            if task.info.status in ["started", "submitted"]:
                execution_time = f"*{execution_time}"
    else:
        execution_time = "n/a"
    return execution_time


def resource_fmtr(executer, unused_task):
    resource_type = f"{executer.host_type} {executer.vm_type}"
    if executer.n_mpi_hosts > 1:
        resource_type += f" x{executer.n_mpi_hosts}"
    return resource_type


class TasksTabulator(BaseTabulator):
    """Tabulator for lists of tasks."""
    default_header_formatter = header_formatter
    id = Col("ID", "id")
    simulator = Col("Simulator", "get_simulator_name")
    status = Col("Status", "info.status", value_formatter=status_formatter)
    submitted = Col("Submitted",
                    "info.input_submit_time",
                    value_formatter=datetime_formatter)
    started = Col("Started",
                  "info.start_time",
                  value_formatter=datetime_formatter)
    computetime = Col("Computation Time", "", value_formatter=computetime_fmtr)
    restype = Col("Resource Type",
                  "info.executer",
                  value_formatter=resource_fmtr)
