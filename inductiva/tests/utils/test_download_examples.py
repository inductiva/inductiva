"""Test examples download."""
import os
import ssl
import sys
import shutil
import zipfile
import platform
import inductiva

_URL = "https://storage.googleapis.com/inductiva-api-demo-files/" \
          "openfoam-input-example.zip"


def test_download_from_url_unzip_false():

    expected_download_path = os.path.join(os.curdir,
                                          "openfoam-input-example.zip")
    # Check that the file does not exist yet.
    assert not os.path.exists(expected_download_path)

    # exposes windows native system certificate stores
    # Solves issues with SSLCertVerificationError
    # truststore is only available for python 3.10 and above
    if platform.system() == "Windows" and sys.version_info >= (3, 10):
        import truststore  # pylint: disable=C0415
        truststore.inject_into_ssl()
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


def test_download_from_url_unzip_true():
    # Disable SSL verification.
    # Had to do this for tests to pass on windows.
    #pylint: disable=protected-access
    ssl._create_default_https_context = ssl._create_unverified_context

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
