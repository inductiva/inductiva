
# Troubleshooting Installation Issues

Encountered unexpected problems while trying to install something? This page
provides a list of common issues and their solutions.

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

---

With these steps, you should be able to resolve the issues and continue with the installation process.
