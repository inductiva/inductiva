"""
Sample usage of the inductiva package.
"""
import inductiva

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.update_config(address="http://192.168.1.50:8000", output_dir="output")

    result = inductiva.math.sum(a=1, b=1)

    logging.debug("Result = %s", result)
