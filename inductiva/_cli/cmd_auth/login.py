"""Register CLI command for login."""
import argparse
import getpass
import os

import inductiva
from inductiva import constants, users, utils

from jose import jwe
from jose.exceptions import JWEError


def is_valid_jwe(token):
    try:
        # Attempt to decode the JWE token
        jwe.get_unverified_header(token)
        return True
    except JWEError:
        return False


def login(args):
    """
    Prompts the user to enter their API key and stores it securely.

    The function will prompt the user to enter their API key, which can be
    obtained from their account at https://console.inductiva.ai/account.
    """

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
            "    Setting a new API key will erase the existing one.")

    prompt_func = getpass.getpass if args.private else input
    warning = " (input will not be visible)" if args.private else ""

    prompt = prompt_func(
        "    To log in, you need an API key. You can obtain it "
        "from your account at https://console.inductiva.ai/account.\n"
        f"Please paste your API key here{warning}: ")

    api_key = prompt.strip()

    if not api_key:
        print("Error: API key cannot be empty.")
        return

    if not is_valid_jwe(api_key):
        print("Error: Invalid API key format.")
        return

    # Set the API key to check if it is valid
    inductiva.set_api_key(api_key)

    # If the API key is invalid, this will raise an exception
    user_info = users.get_info()

    utils.set_stored_api_key(api_key)

    user_name = user_info["name"] or ""
    print(f"Welcome back {user_name}!")


def register(parser):
    """Register the login command."""
    subparser = parser.add_parser("login",
                                  help="Login using Inductiva API key.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva login` command allows you to log in using your key.\n"
        "You can obtain your API key from your account at "
        "https://console.inductiva.ai/account.\n")

    subparser.add_argument("--private",
                           action='store_true',
                           help="Hide API Key.")

    subparser.set_defaults(func=login)
