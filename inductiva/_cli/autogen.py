"""
Utility functions for automatic generation and customization
of CLI documentation using sphinx-argparse-cli.
"""

import argparse
from inductiva._cli.main import get_main_parser

def get_subparser(
    parser: argparse.ArgumentParser,
    name: str
) -> argparse.ArgumentParser:
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    subparser = subparsers.choices[name]
    subparsers.choices.clear()
    return subparser

def get_parser(command: str):
    sub_commands = command.split(" ")
    parser = get_main_parser()
    for sub_command in sub_commands:
        parser = get_subparser(parser, sub_command)
    return parser

def auth_login_parser():
    return get_parser("auth login")

def auth_logout_parser():
    return get_parser("auth logout")