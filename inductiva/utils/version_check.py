class VersionCheckException(Exception):

    def __init__(self, message: str):
        super().__init__(message)


def compare_versions(inductiva_version: str, backend_version: str):
    """
    Compares the installed client version (`inductiva_version`) with the backend version (`backend_version`).

    The function compares the major, minor, and patch versions of both input versions, following the semantic versioning convention.
    If the client version is older than the backend version, it raises a `VersionCheckException`.
    If the client's minor and/or patch versions are behind the backend's, but the major version matches, 
    it prints a message prompting the user to update, but does not raise an exception.

    Args:
        inductiva_version (str): Installed client version to be compared, in the format 'X.Y.Z'.
        backend_version (str): Backend's version to be compared against, in the format 'X.Y.Z'.

    Raises:
        VersionCheckException: If `inductiva_version` is older than `backend_version` in terms of major version.

    Examples:
        >>> compare_versions("1.0.0", "2.0.0")
        VersionCheckException: The installed version (1.0.0) is too old. Please upgrade to version 2.0.0 or newer.

        >>> compare_versions("1.0.0", "1.1.0")
        Update available 1.0.0 -> 1.1.0
        To update to the latest version using pip, run:
        pip install --upgrade inductiva
    """
    inductiva_major, inductiva_minor, inductiva_patch = map(
        int, inductiva_version.split('.'))

    backend_major, backend_minor, backend_patch = map(
        int, backend_version.split('.'))

    if (inductiva_major < backend_major):
        raise VersionCheckException(
            f"The installed version ({inductiva_version}) is too old. Please upgrade to version {backend_version} or newer."
        )

    if inductiva_major != backend_major:
        return

    # Compare minor and patch versions
    if inductiva_minor < backend_minor or (inductiva_minor == backend_minor and
                                           inductiva_patch < backend_patch):

        print(f"Update available {inductiva_version} -> {backend_version}")
        print("To update to the latest version using pip, run:")
        print("pip install --upgrade inductiva")
