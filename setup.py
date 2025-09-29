"""Setup file."""

import os
from setuptools import setup
from setuptools.config import read_configuration

# Read configuration from setup.cfg
config = read_configuration("setup.cfg")
base_install_requires = config["options"]["install_requires"]

setup(install_requires=base_install_requires +
      ([] if os.getenv("NOAIORTC", False) else ["aiortc"]),)
