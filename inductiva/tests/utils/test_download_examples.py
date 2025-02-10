"""Test examples download."""
import os
import ssl
import zipfile
import inductiva


def test_download_from_url():
    # Disable SSL verification.
    # Had to do this for tests to pass on windows.
    #pylint: disable=protected-access
    ssl._create_default_https_context = ssl._create_unverified_context

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
