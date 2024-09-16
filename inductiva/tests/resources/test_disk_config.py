"""Tests for DiskConfig Class"""
import inductiva

import pytest


@pytest.mark.parametrize(
    "max_size_gb, is_resizable",
    [
        (1.5, True),
        (0, True),
        (-1, True),
        (100, "False"),
    ],
)
def test_disk_config_exceptions(max_size_gb, is_resizable):
    with pytest.raises(ValueError) as _:
        _ = inductiva.resources.DiskConfig(max_size_gb=max_size_gb,
                                           is_resizable=is_resizable)


def test_disk_config_not_resizable():
    disk_config = inductiva.resources.DiskConfig(max_size_gb=100,
                                                 is_resizable=False)
    assert disk_config.max_size_gb == 100
    assert not disk_config.is_resizable
    assert isinstance(disk_config.resize_trigger_gb,
                      int) and disk_config.resize_trigger_gb > 0
    assert isinstance(disk_config.resize_increment_gb,
                      int) and disk_config.resize_increment_gb > 0


def test_disk_config_resizable():
    disk_config = inductiva.resources.DiskConfig(max_size_gb=100,
                                                 is_resizable=True)
    assert disk_config.max_size_gb == 100
    assert disk_config.is_resizable
    assert isinstance(disk_config.resize_trigger_gb,
                      int) and disk_config.resize_trigger_gb > 0
    assert isinstance(disk_config.resize_increment_gb,
                      int) and disk_config.resize_increment_gb > 0
