# Install the Inductiva Python Package

Set up in seconds – Install the Inductiva package and start running simulations effortlessly!

### Before You Start

Make sure you have a compatible version of Python (>=3.9) and pip installed, plus some basic Python knowledge.

Check our <a href="/guides/systemrequirements">System Prep guide</a> to help you check these essentials before getting started.

<!-- Check our <a href="https://docs.inductiva.ai/en/latest/preinstallation/system/system-requirements.html#">System Prep guide</a> to help you check these essentials before getting started.   -->

## Step 1: Install the Package

Open your Terminal (Linux/MacOS) or Command Prompt/PowerShell (Windows) and enter:

```python
pip install inductiva
```

## Step 2: Authenticate Using Inductiva's CLI (Command Line Interface)

Now that the Inductiva package is installed, run the authentication command:

```python
inductiva auth login
```

> [!WARNING]
> Windows users might experience an error when trying to authenticate this way, meaning that Inductiva's CLI (Command Line Interface) was not successfuly installed.  
> In that case, follow the instructions to Authenticate Using the Python API. Otherwise, move to Step 3.

### Authenticate Using the Python API

If you weren't able to authenticate using Inductiva's CLI (Command Line Interface), you can do it directly within your Python script.  

On your Terminal (Linux/MacOS) or Command Prompt/PowerShell (Windows) start your Python interpreter:

```python
python
```

You should get a message with the Python version.  
Then type:

```python
import inductiva
inductiva.auth.login()
```

## Step 3: Authenticate With Your API Key

Regardless if you used Inductiva's CLI or the Python API, you should now be getting a prompt "Please paste your API Key here:"  

Retrieve your API Key from [Inductiva's web Console](https://console.inductiva.ai/account/details), and paste it in the Terminal.

To confirm your authentication, type:

```python
inductiva user info
```

This will display your account information, confirming that the API key has been stored successfully.
