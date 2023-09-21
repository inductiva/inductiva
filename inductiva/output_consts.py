"""Connect tasks to simulation output classes via method_name."""
from inductiva import fluids, molecules, coastal

OUTPUT_CONSTS = {
    "wind_tunnel": {
        "output_class":
            fluids.WindTunnelOutput,
        "default_files": [
            "pressure_field.vtk", "streamlines.vtk", "stdout.txt", "stderr.txt",
            "force_coefficients.csv", "xy_flow_slice.vtk", "xz_flow_slice.vtk",
            "yz_flow_slice.vtk", "constant/triSurface/object.obj"
        ]
    },
    "wind_terrain": {
        "output_class": fluids.SteadyStateOutput,
        "default_files": None
    },
    "heat_sink": {
        "output_class": fluids.HeatSinkOutput,
        "default_files": None
    },
    "dam_break": {
        "output_class": fluids.SPHSimulationOutput,
        "default_files": None
    },
    "fluid_block": {
        "output_class": fluids.SPHSimulationOutput,
        "default_files": None
    },
    "fluid_tank": {
        "output_class": fluids.FluidTankOutput,
        "default_files": None
    },
    "coastal_area": {
        "output_class": coastal.CoastalAreaOutput,
        "default_files": None
    },
    "protein_solvation": {
        "output_class":
            molecules.ProteinSolvationOutput,
        "default_files": [
            "protein.gro", "compressed_trajectory.xtc", "solvated_protein.tpr",
            "topol.top", "rmsf_values.npy", "stdout.txt", "stderr.txt"
        ]
    },
    "mdwater_box": {
        "output_class":
            molecules.MDWaterBoxOutput,
        "default_files": [
            "trajectory.xtc", "eql.tpr", "stdout.txt", "stderr.txt"
        ]
    },
    "stellarators": {
        "output_class": None,
        "default_files": None
    },
    "deformable_plate": {
        "output_class": None,
        "default_files": None
    }
}
