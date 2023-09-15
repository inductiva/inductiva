import pathlib
import inductiva


def test_download_from_rcsb(tmp_path: pathlib.Path):
    """Check if files with extensions are correctly created."""
    inductiva.working_dir = tmp_path
    inductiva.molecules.examples.download_from_rcsb("1A3N")
    assert (tmp_path / "1A3N.pdb").exists()
