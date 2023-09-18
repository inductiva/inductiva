"""Constants for connecting tasks to their simulation output classes."""
import inductiva

SCENARIO_METHODS_NAMES = [
    "wind_tunnel", "wind_terrain", "heat_sink", "dam_break", "fluid_block",
    "fluid_tank", "coastal_area", "protein_solvation", "mdwater_box",
    "deformable_plate", "stellarators"
]

OUTPUT_CLASSES = {
    "wind_tunnel": inductiva.fluids.WindTunnelOutput,
    "wind_terrain": inductiva.fluids.SteadyStateOutput,
    "heat_sink": inductiva.fluids.HeatSinkOutput,
    "dam_break": inductiva.fluids.SPHSimulationOutput,
    "fluid_block": inductiva.fluids.SPHSimulationOutput,
    "fluid_tank": inductiva.fluids.FluidTankOutput,
    "coastal_area": inductiva.coastal.CoastalAreaOutput,
    "protein_solvation": inductiva.molecules.ProteinSolvationOutput,
    "mdwater_box": inductiva.molecules.MDWaterBoxOutput,
    "stellarators": None,
    "deformable_plate": None
}

DEFAULT_OUTPUT_FILES = {
    "wind_tunnel": [
        "pressure_field.vtk", "streamlines.vtk", "stdout.txt", "stderr.txt",
        "force_coefficients.csv", "xy_flow_slice.vtk", "xz_flow_slice.vtk",
        "yz_flow_slice.vtk", "constant/triSurface/object.obj"
    ],
    "wind_terrain": None,
    "heat_sink": None,
    "dam_break": None,
    "fluid_block": None,
    "fluid_tank": None,
    "coastal_area": None,
    "protein_solvation": [
        "protein.gro", "compressed_trajectory.xtc", "solvated_protein.tpr",
        "topol.top", "rmsf_values.npy"
    ],
    "mdwater_box": None,
    "stellarators": None,
    "deformable_plate": None
}
