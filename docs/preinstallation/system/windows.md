# Windows System Prep Guide

## TL;DR
Pre-Install Must-Dos:
- [Step 1: Check if Python is Installed](#step-1-check-if-python-is-installed)
- [Step 2: Update pip and Set Up Your PATH for Python](#step-2-update-pip-and-set-up-your-path-for-python)


## Step 1: Check if Python is Installed

First things first, let’s make sure **Python 3** is set up correctly on your system. 
If you already have it installed, great! We’ll be able to move straight to installing 
the Inductiva API like any other Python package.

If it’s missing, we’ll go over how to install it on Windows.

**Step-by-Step**

1. **Open the Command Prompt** 
	
    Click on the Windows icon in the bottom-left corner, type `Command Prompt`, and press **Enter** to open the Command Prompt app.
     
2. **Check for Python 3**
    
    In the Command Prompt, type:
   
    ```bash
    python3
    ```

- *If Python 3 is installed*, you’ll see something like this:

    ```bash
    Python 3.12.7 (tags/v3.12.7:abcdef, Oct  3 2023, 12:00:00) [MSC v.1928 64 bit (AMD64)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    ```

- *If Python 3 isn’t installed*, running python or python3 in the Command Prompt will usually display a message indicating it’s not recognized as an internal or external command.

    Since recent versions of Windows (Windows 10 and later) include a feature to assist with installing Python, you may see something like this:

    ```bash
    'python' is not recognized as an internal or external command, operable program, or batch file.
    ```

    Windows will then prompt you with an option to Open the Microsoft Store, where it takes you directly to the Python page. From there, you can click **Get** to install Python. 

````{eval-rst}
.. important::
   During installation, make sure to check the box that says Add Python.exe to PATH. This ensures the PATH is set automatically, so you won’t have to do it manually later.
````

The installation can take a few minutes to complete.

4. **Test Again**

    Once installation is complete, type `python3` in the Command Prompt again. You should now see confirmation that Python 3 is ready to go!

     ```bash
    Python 3.12.7 (tags/v3.12.7:abcdef, Oct  3 2023, 12:00:00) [MSC v.1928 64 bit (AMD64)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    ```

## Step 2: Update pip and Set Up Your PATH for Python

Now that you have **Python 3** installed, let’s make sure pip, the Python package installer, is up to date. This will help avoid any compatibility issues.

### Update pip

To update pip, open your Command Prompt and type:

```bash
python3 -m pip install --upgrade pip
```

Since we installed Python from the Microsoft Store and selected the Add Python.exe to PATH option, this command should run smoothly without any issues.

---

Awesome! Now that pip is updated and Python is set in your PATH, you’re all set for [installing the Inductiva Python Package](https://console.inductiva.ai/) and start simulating!

If you run into any issues or challenges while installing the API, please reach out to us at support@inductiva.ai. We’d love to help troubleshoot and find ways to make the setup process even smoother.

You can also check out our [troubleshooting guide](../api_reference/troubleshooting.md) for more information.
