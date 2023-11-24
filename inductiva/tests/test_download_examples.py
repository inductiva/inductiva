"""Test examples download."""
import os
import zipfile
import inductiva


def test_download_from_rcsb():
    """Check if files with extensions are correctly created."""

    file_path = inductiva.molecules.utils.download_pdb_from_rcsb("1A3N")
    pdb_file = file_path.rsplit("/", maxsplit=1)[-1]
    assert os.path.exists(file_path)
    assert pdb_file == "1A3N.pdb"


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
