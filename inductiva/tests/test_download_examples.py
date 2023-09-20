"""Test examples download."""
import pathlib
import inductiva


def test_download_from_rcsb(tmp_path: pathlib.Path):
    """Check if files with extensions are correctly created."""
    inductiva.working_dir = tmp_path
    inductiva.molecules.utils.download_pdb_from_rcsb("1A3N")
    assert (tmp_path / "1A3N.pdb").exists()


def test_download_zip(tmp_path: pathlib.Path):
    inductiva.working_dir = tmp_path

    downloaded_to = inductiva.utils.files.download_from_url(
        "https://downloads.emodnet-bathymetry.eu/high_resolution/"
        "590_HR_Lidar_Algarve.emo.zip")

    path = tmp_path / "590_HR_Lidar_Algarve.emo"
    assert str(path.absolute()) == downloaded_to
    assert path.exists()
