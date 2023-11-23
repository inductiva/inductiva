"""Test examples download."""
import os
import pathlib
import zipfile
import inductiva


def test_download_from_rcsb(tmp_path: pathlib.Path):
    """Check if files with extensions are correctly created."""
    inductiva.working_dir = tmp_path
    inductiva.molecules.utils.download_pdb_from_rcsb("1A3N")
    assert (tmp_path / "1A3N.pdb").exists()


def test_download_from_url():
    url = "https://storage.googleapis.com/inductiva-api-demo-files/" \
          "openfoam-input-example.zip"

    # Check if the unzip argument works.
    download_path = inductiva.utils.files.download_from_url(url, unzip=True)
    assert os.path.exists(download_path)
    assert not zipfile.is_zipfile(download_path)

    # Check if the default is a zip file.
    download_path = inductiva.utils.files.download_from_url(url)
    assert os.path.exists(download_path)
    assert zipfile.is_zipfile(download_path)
