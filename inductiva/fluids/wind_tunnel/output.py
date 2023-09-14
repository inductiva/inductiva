"""Visualization processing of WindTunnel scenario.
Currently, we only support the OpenFOAM simulator.
"""
import os
import csv

from inductiva import fluids, types


class WindTunnelOutput(fluids.SteadyStateOutput):
    """Post-Process WindTunnel simulation outputs.
    This class is based on a more general class for post-processing the
    outputs of steady-state simulations. It inherits the following tools:
        - Pressure over object;
        - Cutting plane;
        - Stream lines;

    Additionally, it provides the following methods:
        - Force coefficients (added here only for this scenario).

    Current Support:
        OpenFOAM.
    """

    def get_force_coefficients(self, save_path: types.Path = None):
        """Get the force coefficients of the object in the WindTunnel.
        
        The force coefficients are provided in a .dat file during the
        simulation run-time. This file contains 8 lines that are provide
        the general input information. In this function, we read the file,
        ignore the first 8 lines and read the force coefficients for the 
        simulation_time chosen.

        Args:
            save_path: Path to save the force coefficients in a .csv file.
        """

        if not self.full_output:
            force_coefficients_path = os.path.join(self.sim_output_path,
                                                   "force_coefficients.csv")
            with open(force_coefficients_path, "r",
                      encoding="utf-8") as csv_file:
                csv_reader = csv.reader(csv_file)
                force_coefficients = list(csv_reader)
        else:
            num_header_lines = 8
            force_coefficients_path = os.path.join(self.sim_output_path,
                                                   "postProcessing",
                                                   "forceCoeffs1", "0",
                                                   "forceCoeffs.dat")
            force_coefficients = []
            last_iteration = int(self.get_last_iteration())

            with open(force_coefficients_path, "r",
                      encoding="utf-8") as forces_file:
                for index, line in enumerate(forces_file.readlines()):
                    # Pick the line 8 of the file:
                    # [#, Time, Cm, Cd, Cl, Cl(f), Cl(r)] and remove # column
                    if index == num_header_lines:
                        force_coefficients.append(line.split()[1:])
                    # Add the force coefficients for the simulation time chosen
                    elif index == num_header_lines + last_iteration + 1:
                        force_coefficients.append(line.split())

            if save_path:
                with open(save_path, "w", encoding="utf-8") as csv_file:
                    csv_writer = csv.writer(csv_file)
                    csv_writer.writerows(force_coefficients)

        return force_coefficients
