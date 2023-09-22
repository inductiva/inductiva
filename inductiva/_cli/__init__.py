"""Inductiva CLI.

`main` is the entrypoint for the CLI. It creates a parser and registers
subparsers for several subcommands (e.g. `tasks` and `resources`).
"""
from . import utils
from .main import main
from .tasks import register_tasks_cli
from .resources import register_resources_cli
