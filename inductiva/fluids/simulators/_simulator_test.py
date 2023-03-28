"""Test Simulator class."""
import pathlib
import pytest
from inductiva.fluids.simulators.swash import SWASH


def test_simulator_constructor(tmp_path: pathlib.Path):
    """Test constructor of the Simulator class."""

    config_filename = "config.json"
    config_file = tmp_path / config_filename
    config_file.write_text("{}")

    invalid_dir = tmp_path / "does_not_exist"

    with pytest.raises(ValueError):
        SWASH(invalid_dir, "sim_config_filename")

    with pytest.raises(ValueError):
        SWASH(tmp_path, "filename")

    SWASH(tmp_path, config_filename)
