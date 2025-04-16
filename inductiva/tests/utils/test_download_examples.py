"""Test examples download."""
import os
import shutil
import ssl
import urllib.request
import zipfile
from unittest import mock

import inductiva

_URL = "https://storage.googleapis.com/inductiva-api-demo-files/" \
          "openfoam-input-example.zip"


@mock.patch("urllib.request.urlopen",
            side_effect=lambda *args, **kwargs: urllib.request.urlopen(
                args[0], context=ssl._create_unverified_context()))
def test_download_from_url_unzip_false(mock_urlopen):

    expected_download_path = os.path.join(os.curdir,
                                          "openfoam-input-example.zip")
    # Check that the file does not exist yet.
    assert not os.path.exists(expected_download_path)
    # Download file from url.
    download_path = inductiva.utils.files.download_from_url(_URL)
    # Check if the file exists.
    assert os.path.exists(download_path)
    # Check if the file was downloaded to the expected location.
    assert os.path.samefile(download_path, expected_download_path)
    # Check if the file is a zip file.
    assert zipfile.is_zipfile(download_path)
    # Delete the zip file that was downloaded.
    os.remove(download_path)


@mock.patch("urllib.request.urlopen",
            side_effect=lambda *args, **kwargs: urllib.request.urlopen(
                args[0], context=ssl._create_unverified_context()))
def test_download_from_url_unzip_true(mock_urlopen):
    # Disable SSL verification.
    # Had to do this for tests to pass on windows.
    #pylint: disable=protected-access
    # ssl._create_default_https_context = ssl._create_unverified_context

    # Check that the expected folder to unzip to, does not exist yet.
    expected_unzipped_folder_path = os.path.join(os.curdir,
                                                 "openfoam-input-example")
    assert not os.path.exists(expected_unzipped_folder_path)

    # Download and decompress file from url.
    download_path = inductiva.utils.files.download_from_url(_URL, unzip=True)

    # Check if the unzipped folder has now been created.
    assert os.path.isdir(download_path)

    # Check folder was decompressed to the expected location.
    assert os.path.samefile(download_path, expected_unzipped_folder_path)

    # Delete the unzipped folder.
    shutil.rmtree(download_path)
