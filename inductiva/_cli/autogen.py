"""
Utility functions for automatic generation and customization
of CLI documentation using sphinx-argparse-cli.
"""

import argparse
from typing import Optional
from inductiva._cli.main import get_main_parser

def get_subparsers(
    parser: argparse.ArgumentParser,
) -> Optional[argparse._SubParsersAction]:
    return next(
        (action for action in parser._actions
         if isinstance(action, argparse._SubParsersAction)),
        None,
    )

def get_subparser(
    parser: argparse.ArgumentParser,
    name: str
) -> argparse.ArgumentParser:
    subparsers = get_subparsers(parser)
    return subparsers.choices[name]

def get_parser(command: str):
    sub_commands = command.split(" ")
    parser = get_main_parser()
    for sub_command in sub_commands:
        parser = get_subparser(parser, sub_command)
    subparsers = get_subparsers(parser)
    if subparsers:
        subparsers.choices.clear()
    return parser

def auth_parser():
    return get_parser("auth")

def auth_login_parser():
    return get_parser("auth login")

def auth_logout_parser():
    return get_parser("auth logout")