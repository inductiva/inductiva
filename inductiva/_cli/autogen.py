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
    return subparsers.choices[name]

def auth_parser():
    main_parser = get_main_parser()
    subparser = get_subparser(main_parser, "auth")
    subparsers = next(
        action for action in subparser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    subparsers.choices.clear()
    return subparser


def auth_login_parser():
    main_parser = get_main_parser()
    auth_parser = get_subparser(main_parser, "auth")
    return get_subparser(auth_parser, "login")

def auth_logout_parser():
    main_parser = get_main_parser()
    auth_parser = get_subparser(main_parser, "auth")
    return get_subparser(auth_parser, "logout")
