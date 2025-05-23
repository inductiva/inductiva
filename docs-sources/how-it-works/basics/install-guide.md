# Install the Inductiva Python Package

Set up in seconds – Install the Inductiva package and start running simulations effortlessly!

### Before You Start

Make sure you have a compatible version of Python (>=3.9) and pip installed, plus some basic Python knowledge.

Check our [System Prep guide](https://docs.inductiva.ai/en/latest/preinstallation/system/system-requirements.html#) to help you check these essentials before getting started.

<!-- Check our <a href="https://docs.inductiva.ai/en/latest/preinstallation/system/system-requirements.html#">System Prep guide</a> to help you check these essentials before getting started.   -->

## Step 1: Install the Package

Open your Terminal (Linux/MacOS) or Command Prompt/PowerShell (Windows) and enter:

```python
pip install inductiva
```

## Step 2: Authenticate With Your API Key

You have multiple ways to authenticate with the Inductiva Python package:

### Option 1: Using the Command Line

In your Command Prompt, run the following authentication command:

```python
inductiva auth login
```

When prompted, paste your unique API key, which you can retrieve from [Inductiva's web Console](https://console.inductiva.ai/account/details).

To confirm your authentication and view user details, run:

```python
inductiva user info
```

This will display your account information, confirming that the API key has been stored successfully.


### Option 2: Using the Python API

Alternatively, you can authenticate directly within your Python script:

```python
import inductiva

inductiva.auth.login()
```
When prompted, paste your unique API key, which you can retrieve from [Inductiva's web Console](https://console.inductiva.ai/account/details).
