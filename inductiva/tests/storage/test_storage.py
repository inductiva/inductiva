import pytest
from unittest import mock
from inductiva.storage.storage import _construct_remote_paths


@pytest.mark.parametrize(
    "local_path, remote_dir, is_dir, mock_list_files_return, expected_remote_paths",
    [
        # Case 1: Directory (Unix style)
        (
            '/local/path/to/dir',
            'remote/dir',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 2: File (Unix style)
        (
            '/local/path/to/file.txt',
            'remote/dir',
            False,
            ([], 0),
            ['remote/dir/file.txt'],
        ),
        # Case 3: Trailing slash in remote_dir (Unix)
        (
            '/local/path/to/dir',
            'remote/dir/',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 4: Directory (Windows style)
        (
            r'C:\\local\\path\\to\\dir',
            r'remote\\dir',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 5: File (Windows style)
        (
            r'C:\\local\\path\\to\\file.txt',
            r'remote\\dir',
            False,
            ([], 0),
            ['remote/dir/file.txt'],
        ),
        # Case 6: Trailing slash in remote_dir (Windows style)
        (
            r'C:\\local\\path\\to\\dir',
            r'remote\\dir\\',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 7: Trailing slash in remote_dir (Windows style)
        (
            r'local\\path\\to\\dir',
            r'remote\\dir\\',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 8: Relative directory path (Windows style)
        (
            r'\\relative\\path\\to\\dir',
            r'remote\\dir',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 9: Relative file path (Windows style)
        (
            r'\\another\\path\\file.txt',
            r'remote\\dir',
            False,
            ([], 0),
            ['remote/dir/file.txt'],
        )
    ])
@mock.patch('inductiva.storage.storage._list_files')
@mock.patch('os.path.isdir')
@mock.patch('os.path.getsize', return_value=0)
def test__construct_remote_paths(
    _,
    mock_isdir_func,
    mock_list_files_func,
    local_path,
    remote_dir,
    is_dir,
    mock_list_files_return,
    expected_remote_paths,
):

    mock_isdir_func.return_value = is_dir
    mock_list_files_func.return_value = mock_list_files_return

    remote_file_paths, _, _ = _construct_remote_paths(
        local_path,
        remote_dir,
    )

    assert remote_file_paths == expected_remote_paths
