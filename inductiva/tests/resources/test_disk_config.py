import inductiva


def test_disk_config_only_size():
    disk = inductiva.resources.DiskConfig(size_gb=1)
    assert disk.size_gb == 1


def test_disk_config_all_arguments():
    disk = inductiva.resources.DiskConfig(size_gb=1,
                                          resize_trigger_gb=2,
                                          resize_increment_gb=3,
                                          max_size_gb=4)
    assert (disk.size_gb == 1 and disk.resize_trigger_gb == 2 and
            disk.resize_increment_gb == 3 and disk.max_size_gb == 4)


def test_disk_config_resize_config():
    disk = inductiva.resources.DiskConfig(size_gb=1,
                                          resize_trigger_gb=2,
                                          resize_increment_gb=3,
                                          max_size_gb=4)
    assert disk.resize_config == {
        "free_space_threshold_gb": 2,
        "size_increment_gb": 3,
        "max_disk_size_gb": 4,
    }


def test_disk_config_resize_config_not_resizable():
    disk = inductiva.resources.DiskConfig(size_gb=1)
    assert disk.resize_config is None


def test_disk_config_not_resizable():
    disk = inductiva.resources.DiskConfig(size_gb=1)
    assert disk.is_resizable is False


def test_disk_config_is_resizable():
    disk = inductiva.resources.DiskConfig(size_gb=1,
                                          resize_trigger_gb=2,
                                          resize_increment_gb=3,
                                          max_size_gb=4)
    assert disk.is_resizable is True
