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
    "inductiva_sph @ git+https://github.com/inductiva/inductiva-sph.git@673f06311dd78846c7c2ed3b7344e292f85868eb",
    "inductiva_data @ git+https://github.com/inductiva/inductiva-data.git@40b32c0207a828ba976eed1ad490984c1985e911",
    "inductiva_utils @ git+https://github.com/inductiva/inductiva-utils.git",
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
