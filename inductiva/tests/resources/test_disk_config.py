"""Tests for DiskConfig Class"""
import inductiva


def test_disk_config_only_size():
    disk = inductiva.resources.DiskConfig(max_size_gb=100, is_resizable=False)
    assert disk.max_size_gb == 100
    assert disk.is_resizable is False
    assert disk.resize_config is None


def test_disk_config_resizable():
    disk = inductiva.resources.DiskConfig(max_size_gb=100)
    assert disk.is_resizable is True

    assert isinstance(disk.resize_config["free_space_threshold_gb"], int)
    assert isinstance(disk.resize_config["size_increment_gb"], int)
    assert disk.resize_config["max_disk_size_gb"] == 100
    assert isinstance(disk.resize_config, dict)


def test_disk_config_resizable_all_arguments():
    disk = inductiva.resources.DiskConfig(max_size_gb=100, is_resizable=True)
    assert disk.is_resizable is True

    assert isinstance(disk.resize_config["free_space_threshold_gb"], int)
    assert isinstance(disk.resize_config["size_increment_gb"], int)
    assert disk.resize_config["max_disk_size_gb"] == 100
    assert isinstance(disk.resize_config, dict)
