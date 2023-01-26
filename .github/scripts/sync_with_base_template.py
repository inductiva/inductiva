"""Script for syncing the current repository with `inductiva/start-a-project`.

DANGER!!! THIS SCRIPT MAY DELETE YOUR LOCAL CHANGES. USE ON CLEAN BRANCHES ONLY.

Many projects started from `inductiva/start-a-project` directly; others were
later migrated to follow the same project structure. Over time, though, we
realized we wanted to make more changes to `inductiva/start-a-project` and have
an easy way to propagate them to the other repositories. This script tries to
accomplish exactly that.

After synchronization, changes have to be inspected and committed manually.
"""

import os
import sys
import tempfile

FILES = [
    ".gitattributes",
    ".github/workflows",
    ".gitignore",
    ".pylintrc",
    ".readthedocs.yml",
    ".style.yapf",
    "conda.yml",
    "docs/conf.py",
    "docs/requirements.txt",
]

if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="inductiva-") as tmpdir:
        # Clone base repository.
        os.system(
            f"git clone https://github.com/inductiva/start-a-project {tmpdir}")

        # Get latest list of files to sync.
        sys.path.insert(0, os.path.join(tmpdir, ".github/scripts/"))
        from sync_with_base_template import FILES as FILES_TO_SYNC

        # Sync files.
        for file in FILES_TO_SYNC:
            src = os.path.join(tmpdir, file)
            subdir = os.path.dirname(file)
            # creates subdirectory (e.g. docs/) as it might not exist
            if subdir and not os.path.exists(subdir):
                os.system(f"mkdir -p '{subdir}'")
            os.system(f"rm -rf '{file}'")
            os.system(f"cp -R '{src}' '{file}'")

        # Sync current file.
        file = ".github/scripts/sync_with_base_template.py"
        src = os.path.join(tmpdir, file)
        os.system(f"rm -f '{file}'")
        os.system(f"cp '{src}' '{file}'")
