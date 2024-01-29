"""Inductiva CLI.

`main` is the entrypoint for the CLI. It creates a parser and registers
subparsers for several subcommands (e.g. `tasks` and `machines`).
"""
from . import utils
from . import machines
from .main import main
from .tasks import register_tasks_cli
from .machines_parsers import register_machines_cli
from .logs import register_logs_cli
