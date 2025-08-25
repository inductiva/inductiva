""" Unit tests for the storage module."""

import pytest
from unittest import mock
from inductiva import utils
from inductiva.storage.storage import _construct_remote_paths
from inductiva.utils.data import _normalize_file


@pytest.mark.parametrize(
    'local_path, remote_dir, is_dir, mock_list_files_return,' \
    'expected_remote_paths',
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
            'C:\\local\\path\\to\\dir',
            'remote\\dir',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 5: File (Windows style)
        (
            'C:\\local\\path\\to\\file.txt',
            'remote\\dir',
            False,
            ([], 0),
            ['remote/dir/file.txt'],
        ),
        # Case 6: Trailing slash in remote_dir (Windows style)
        (
            'C:\\local\\path\\to\\dir',
            'remote\\dir\\',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 7: Trailing slash in remote_dir (Windows style)
        (
            'local\\path\\to\\dir',
            'remote\\dir\\',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 8: Relative directory path (Windows style)
        (
            '\\relative\\path\\to\\dir',
            'remote\\dir',
            True,
            (['file1.txt', 'file2.txt'], 100),
            ['remote/dir/file1.txt', 'remote/dir/file2.txt'],
        ),
        # Case 9: Relative file path (Windows style)
        (
            '\\another\\path\\file.txt',
            'remote\\dir',
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


@pytest.mark.parametrize('url,expected', [
    ('http://example.com/my%20file.txt', '/my file.txt'),
    ('https://server/folder/path%3Cname%3E.zip', '/folder/path<name>.zip'),
])
def test_unquote_url_path(url, expected):
    assert utils.data.unquote_url_path(url) == expected


@pytest.mark.parametrize('raw,expected', [
    ('un:safe<name>|file?.txt', 'un_safe_name__file_.txt'),
    ('already_safe.txt', 'already_safe.txt'),
])
def test_sanitize_path(raw, expected):
    assert utils.data.sanitize_path(raw) == expected


@pytest.mark.parametrize(
    'content, suffix, expected',
    [
        # Case 1: CRLF endings in .txt → should normalize to LF
        (b'line1\r\nline2\r\nline3\r\n', '.txt', b'line1\nline2\nline3\n'),
        # Case 2: LF endings in .txt → unchanged
        (b'line1\nline2\nline3\n', '.txt', b'line1\nline2\nline3\n'),
        # Case 3: Non .txt/.sh file → unchanged
        (b'line1\r\nline2\r\n', '.bin', b'line1\r\nline2\r\n'),
        # Case 4: CRLF endings in .sh → should normalize to LF
        (b'#!/bin/bash\r\necho Hello\r\n', '.sh', b'#!/bin/bash\necho Hello\n'),
    ])
def test_normalize_file(tmp_path, content, suffix, expected):
    tmp_file = tmp_path / f'test{suffix}'
    tmp_file.write_bytes(content)

    _normalize_file(str(tmp_file))

    result = tmp_file.read_bytes()

    assert result == expected
