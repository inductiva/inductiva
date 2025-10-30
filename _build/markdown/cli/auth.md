# auth

## inductiva auth - CLI interface

```default
inductiva auth [-h] {login,logout} ...
```

Authentication management utilities.

The `inductiva auth` command allows you to manage your authentication on Inductiva API.

This is a fundamental aspect of your experience: you will only be able  to start machines and launch simulations after you are authenticated.  Once authenticated, your credentials will be stored locally for future  sessions.

However, you will need to authenticate separately on each local machine from which you want to use Inductiva.

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva auth login

```default
inductiva auth login [-h] [--private] [--no-examples]
```

The `inductiva auth login` command logs you in using your API key.

You will be prompted to enter your API key, which you can find in your account on the Web Console at:
    https://console.inductiva.ai/account

Once authenticated, your credentials will be securely stored locally, so you will not need to log in again for future sessions.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--private`**]() - Hide API Key.
* [**`--no-examples`**]() - Does not download example scripts.

#### Examples

```bash

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
```

---

### inductiva auth logout

```default
inductiva auth logout [-h]
```

The `inductiva auth logout` command allows you to log out by removing the stored API key.

This will remove the locally stored API key, requiring you to log in  again for future sessions.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

#### Examples

```bash

$ inductiva auth logout
Logout successful.
```

---
::docsbannersmall
::
