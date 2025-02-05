"""Utility functions for tasks command."""
from typing import Tuple, Optional

from inductiva import tasks


def validate_task_computation_started(
        task: tasks.Task) -> Tuple[bool, Optional[str]]:
    info = task.get_info()
    if info.is_terminal:
        return (False, f"Task {task.id} has terminated.\n"
                "Access its output using:\n\n"
                f"  inductiva tasks download --id {task.id}")
    if not info.status == "computation-started":
        return (False, f"Task {task.id} has not started yet.\n"
                "Wait for computation to start.")

    return (True, None)
