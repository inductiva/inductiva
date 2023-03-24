"""Setup file"""

from setuptools import setup, find_packages  # noqa: H301

NAME = "inductiva"
VERSION = "0.1"

# pylint: disable=line-too-long

REQUIRES = [
    "setuptools >= 21.0.0",
    "absl-py",
    "numpy",
    "scipy",
    "inductiva_web_api_client @ git+https://github.com/inductiva/inductiva-web-api-client.git",
    "inductiva_data @ git+https://github.com/inductiva/inductiva-data.git",
    "inductiva_utils @ git+https://github.com/inductiva/inductiva-utils.git",
    "jinja2",
    "tqdm",
    "gmsh",
    "vtk",
    "xarray",
]

setup(
    name=NAME,
    version=VERSION,
    description="",
    author="Inductiva Research Labs",
    author_email="contact@inductiva.ai",
    url="",
    keywords=["Inductiva"],
    python_requires=">=3.7",
    install_requires=REQUIRES,
    packages=find_packages(exclude=["test", "tests"]),
    include_package_data=True,
    long_description="",
)
