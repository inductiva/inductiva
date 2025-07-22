import argparse
from .main import get_main_parser

def get_subparser(name: str):
    main_parser = get_main_parser()
    subparsers = next(
        action for action in main_parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices[name]

def auth_parser():
    return get_subparser("auth")
