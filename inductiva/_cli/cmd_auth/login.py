"""Register CLI command for login."""
import argparse
import getpass
import os

import inductiva
from inductiva import constants, users, utils


def login(args):
    """Prompts the user to enter their API Key and stores it securely."""

    # pylint: disable=trailing-whitespace,line-too-long
    inductiva_art = r"""     ___  _   _  ____   _   _   ____  _____  ___ __     __ _    
    |_ _|| \ | ||  _ \ | | | | / ___||_   _||_ _|\ \   / // \   
     | | |  \| || | | || | | || |      | |   | |  \ \ / // _ \  
     | | | |\  || |_| || |_| || |___   | |   | |   \ V // ___ \ 
    |___||_| \_||____/  \___/  \____|  |_|  |___|   \_//_/   \_\
    """
    print(inductiva_art)
    if os.path.exists(constants.API_KEY_FILE_PATH):
        print(
            "    You are already logged in. Run `inductiva auth logout` if you "
            "want to log out. \n"
            "    Setting a new API Key will erase the existing one.")

    prompt_func = getpass.getpass if args.private else input
    warning = " (input will not be visible)" if args.private else ""

    prompt = prompt_func(
        "    To log in, you need an API Key. You can obtain it "
        "from your account at https://console.inductiva.ai/account.\n"
        f"Please paste your API Key here{warning}: ")

    api_key = prompt.strip()

    # Set the API Key
    inductiva.set_api_key(api_key, login_message=False)

    # If the API Key is invalid, this will raise an exception
    user_info = users.get_info()

    utils.set_stored_api_key(api_key)

    user_name = user_info["name"] or ""
    print(f"Welcome back {user_name}!")


def register(parser):
    """Register the login command."""
    subparser = parser.add_parser("login",
                                  help="Login using Inductiva API Key.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva login` command allows you to log in using your key.\n"
        "You can obtain your API Key from your account at "
        "https://console.inductiva.ai/account.\n")

    subparser.add_argument("--private",
                           action="store_true",
                           help="Hide API Key.")

    subparser.set_defaults(func=login)
