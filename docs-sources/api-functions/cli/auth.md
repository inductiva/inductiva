# auth

The `inductiva auth` command allows you to manage your
authentication on Inductiva API.

This is a fundamental aspect of your experience: you will only be
able to start machines and launch simulations after you are authenticated.
Ooce authenticated, your credentials will be stored locally for future
sessions. Howecer, you will need to perform the authentication step from
every machine you want to use Inductiva from.

## Usage

```sh
inductiva auth [-h] {login,logout} ...
```

### Options
- **`-h, --help`** → Show help message and exit.

## Available Subcommands

### `login`
Authenticate using your Inductiva API key.

```sh
inductiva auth login
```

You will be prompted to enter your API key.

```sh
     ___  _   _  ____   _   _   ____  _____  ___ __     __ _
    |_ _|| \ | ||  _ \ | | | | / ___||_   _||_ _|\ \   / // \
     | | |  \| || | | || | | || |      | |   | |  \ \ / // _ \
     | | | |\  || |_| || |_| || |___   | |   | |   \ V // ___ \
    |___||_| \_||____/  \___/  \____|  |_|  |___|   \_//_/   \_\

    You are already logged in. Run `inductiva auth logout` if you want to log out.
    Setting a new API Key will erase the existing one.
    To log in, you need an API Key. You can obtain it from your account at https://console.inductiva.ai/account.
Please paste your API Key here:
```

At this stage, please copy paste your personal API key that is available
from your [User Account] (https://console.inductiva.ai/account/profile)
page in the Web Console.  Once authenticated, your credentials will be
stored locally for future sessions using that machine.

### `logout`
Remove stored authentication credentials and log out.

```sh
inductiva auth logout
```

This will remove the locally stored API key, requiring you
to log in again for future interactions.

## Need Help?
Run the following command for more details:

```sh
inductiva auth --help
```
