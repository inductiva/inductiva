"""
Sample usage of the inductiva package.
"""
import inductiva

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    # inductiva.init(address="http://192.168.1.50:8000", output_dir="output")

    simulation = inductiva.fluids.DamBreak(fluid=inductiva.fluids.WATER,
                                           fluid_dimensions=[0.1, 0.2, 0.2])
    simulation_output = simulation.simulate()
    simulation_output.render()
