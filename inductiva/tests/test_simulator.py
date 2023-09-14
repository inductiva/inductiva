"""Tests for the Simulator class."""
import inductiva


def test_override_api_method_prefix():
    simulator = inductiva.simulators.Openfoam()
    assert simulator.api_method_name == "fvm.openfoam.run_simulation"
    simulator.override_api_method_prefix("windtunnel")
    assert simulator.api_method_name == "windtunnel.openfoam.run_simulation"
