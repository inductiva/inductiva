"""Setup file."""

import os
from setuptools import setup
from setuptools.config import read_configuration

# Read configuration from setup.cfg
config = read_configuration("setup.cfg")
base_install_requires = config["options"]["install_requires"]

# Conditionally add aiortc dependency
# aiortc is needed only for the real-time access to output files of a task,
# during its execution. As some Python environments (e.g. pyodide) do not
# support this library, we introduce an environment variable to prevent its
# unnecessary installation. All other features of the inductiva package will
# still work.
setup(install_requires=base_install_requires +
      ([] if os.getenv("NOAIORTC") else ["aiortc"]),)
