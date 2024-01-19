"""Test examples download."""
import os
import zipfile
import inductiva


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
