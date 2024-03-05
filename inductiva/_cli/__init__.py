"""Inductiva CLI.

`main` is the entrypoint for the CLI. It creates a parser and registers
subparsers for several subcommands (e.g. `tasks` and `machines`).
"""
from . import utils
from .main import main, get_main_parser
