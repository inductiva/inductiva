# Troubleshooting Guide
In this troubleshooting section, we aim to guide you through resolving issues
related to incorrect or incomplete installation of the Inductiva package
or its dependencies. This guide assumes you’ve gone through all steps we've shared
in the [user console](https://console.inductiva.ai/), and have faced
some issues.

If you find bugs, need help, or want to talk to the developers, reach out to us on
support@inductiva.ai

If you find any security issues, please report to security@inductiva.ai

---

## Installation Failures
If you encounter issues installing the Inductiva API package, there are several
steps you can take:

### Ensure you have a working pip
```
pip install --upgrade pip
```
### Use virtualenv or venv to isolate dependencies

If installing the package fails, you can retry it on a new Python virtual environment.
A [virtual environment](https://docs.python.org/3/library/venv.html) allows you to
have a fresh Python environment with isolated dependencies.

In your shell, run:

```
python -m venv <venv>
```

In that command, you should replace `<venv>` with the path (*e.g.*, `.venv`) in
which you would like to create the environment. Then, to activate the environment
(again, correctly replacing `<venv>`), run:

For `bash`/`zsh`:

```
source <venv>/bin/activate
```

For `cmd.exe` (Windows):

```
<venv>\Scripts\activate.bat
```

For `PowerShell` (Windows):
```
<venv>\Scripts\Activate.ps1
```

After activating the virtual environment, you can install the package as described
below:

```
pip install --upgrade inductiva
```

---

## Unable to Create `inductiva.exe`

### Problem
You receive an error message like `unable to create inductiva.exe` when running
`pip install inductiva`. This issue is most common on Windows systems and
typically occurs due to insufficient permissions to create the executable file.

### Solution

1. **Run the Command Prompt as Administrator**:
   - Right-click the Command Prompt application.
   - Select **Run as Administrator** from the context menu.

2. **Retry the Installation**:
   - Once the Command Prompt is running with administrative privileges, rerun
   the following command:
     ```bash
     pip install inductiva
     ```

---

## `inductiva` is Not Recognized

### Problem
When attempting to run `inductiva`, you see the error:
```
'inductiva' is not recognized as an internal or external command.
```
This happens for one of two reasons:
1. The `inductiva.exe` file was not created successfully. (Refer to [Unable to Create `inductiva.exe`](#unable-to-create-inductivaexe)).
2. The `Scripts` folder from your Python installation directory is missing from
your PATH environment variable. You might also see a warning like this during
installation:
   ```
   WARNING: The script inductiva.exe is installed in 'C:\path\to\python\Scripts' which is not on PATH.
     Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
   ```

### Solution

#### On Windows:
1. **Open System Properties**:
   - Press `Windows + R` to open the Run dialog.
   - Type `sysdm.cpl` and click **OK**.

2. **Edit Environment Variables**:
   - Go to the **Advanced** tab and click **Environment Variables**.
   - Under **System variables**, find the `PATH` variable and click **Edit**.
   - Click **New** and add the path to the `Scripts` folder in your Python installation directory.  
     Example: `C:\path\to\python\Scripts`

#### On macOS:
1. **Open Your Shell Profile**:
   - Open Terminal and edit your shell profile (e.g., `~/.zshrc`):
     ```bash
     nano ~/.zshrc
     ```

2. **Add the Scripts Folder to the PATH Variable**:
   - Add the following line to the end of the file:
     ```bash
     export PATH="/path/to/python/Scripts:$PATH"
     ```

3. **Restart Terminal**:
   - Restart your terminal for the changes to take effect.

> **Note:** Replace `/path/to/python/Scripts` with the actual path to the
`Scripts` folder on your system. The actual path can be seen in the warning message.
