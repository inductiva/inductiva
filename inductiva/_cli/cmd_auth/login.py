"""Register CLI command for login."""
from typing import List
from tqdm import tqdm
import argparse
import requests
import getpass
import pkgutil
import os
import textwrap

import inductiva
from inductiva import constants, users, utils

INDUCTIVA_ART = \
r"""     ___  _   _  ____   _   _   ____  _____  ___ __     __ _
    |_ _|| \ | ||  _ \ | | | | / ___||_   _||_ _|\ \   / // \
     | | |  \| || | | || | | || |      | |   | |  \ \ / // _ \
     | | | |\  || |_| || |_| || |___   | |   | |   \ V // ___ \
    |___||_| \_||____/  \___/  \____|  |_|  |___|   \_//_/   \_\
"""


def download_example_scripts() -> List[str]:
    """Downloads example scripts from the Inductiva Github repository.

    The scripts are downloaded to the `inductiva_examples` directory.

    Returns:
        list: A list of the names of the downloaded example scripts.
    """
    modules = [
        name
        for _, name, _ in pkgutil.iter_modules(inductiva.simulators.__path__)
    ]
    #openfoam as a different naming convention
    modules.append("openfoam_esi")
    modules.append("openfoam_foundation")

    os.makedirs("inductiva_examples", exist_ok=True)
    downloaded = []
    for module in tqdm(modules,
                       leave=False,
                       desc="Downloading example scripts",
                       unit="script"):
        url = constants.INDUCTIVA_GIT_EXAMPLES_URL + f"{module}/{module}.py"

        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            downloaded.append(f"{module}.py")
            with open(f"inductiva_examples/{module}.py", "w",
                      encoding="utf-8") as f:
                f.write(response.text)
    return downloaded


def login(args):
    """Prompts the user to enter their API Key and stores it securely."""

    # pylint: disable=trailing-whitespace,line-too-long

    print(INDUCTIVA_ART)
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

    first_log_in = utils.set_stored_api_key(api_key)

    user_name = user_info.name or ""

    if first_log_in:
        print(f"Welcome back {user_name}!")
    else:
        print("\n")
        print(f" ■ Welcome {user_name}!\n"
              "Since this is your first time logging in, we will download "
              "some example scripts for you to get started.\n")

        downloaded = download_example_scripts()

        print("The examples are located in the `inductiva_examples` folder.\n"
              "You can run them using the command "
              "`python inductiva_examples/<example.py>`.\n")
        print("Available examples:")
        for i in range(0, len(downloaded), 3):
            line = [word.ljust(25) for word in downloaded[i:i + 3]]
            print("".join(line))
        print("\nRun your first simulation with "
              "`python inductiva_examples/openfoam_esi.py`\n")


def register(parser):
    """Register the login command."""
    subparser = parser.add_parser("login",
                                  help="Login using Inductiva API Key.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = textwrap.dedent("""\
        The `inductiva auth login` command logs you in using your API key.

        You will be prompted to enter your API key, which you can find in
        your account on the Web Console at:
            https://console.inductiva.ai/account

        Once authenticated, your credentials will be securely stored locally,
        so you will not need to log in again for future sessions.
    """)

    subparser.add_argument("--private",
                           action="store_true",
                           help="Hide API Key.")

    subparser.epilog = textwrap.dedent(r"""
        examples:
            $ inductiva auth login
                 ___  _   _  ____   _   _   ____  _____  ___ __     __ _
                |_ _|| \ | ||  _ \ | | | | / ___||_   _||_ _|\ \   / // \\
                 | | |  \| || | | || | | || |      | |   | |  \ \ / // _ \\
                 | | | |\  || |_| || |_| || |___   | |   | |   \ V // ___ \\
                |___||_| \_||____/  \___/  \____|  |_|  |___|   \_//_/   \_\\

                To log in, you need an API Key. You can obtain it from your 
                account at https://console.inductiva.ai/account.
            Please paste your API Key here: 0123456789
            
            ■ Welcome User!
    """)

    subparser.set_defaults(func=login)
