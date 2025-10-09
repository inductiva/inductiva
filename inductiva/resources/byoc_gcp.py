"""GCP BYOC (Bring Your Own Cloud) utilities."""
import subprocess
import re
import tempfile
import os
import datetime
from typing import Optional, Union

from inductiva import users


def estimate_vcpus_from_machine_type(machine_type):
    """Estimate vCPUs from machine type string."""
    # Extract number from machine type (e.g., "c2d-standard-8" -> 8)
    match = re.search(r"-(\d+)$", machine_type)
    if match:
        return int(match.group(1))
    raise ValueError(
        f"Could not extract vCPU count from machine type: {machine_type}")


def create_gcp_startup_script() -> str:
    """Create the startup script for GCP VM."""
    script_path = os.path.join(os.path.dirname(__file__),
                               "gcp_startup_script.txt")

    with open(script_path, "r", encoding="utf-8") as f:
        return f.read()


def check_gcloud_installed() -> bool:
    """Check if gcloud CLI is installed."""
    try:
        subprocess.run(["gcloud", "--version"],
                       capture_output=True,
                       text=True,
                       check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def check_gcloud_auth() -> bool:
    """Check if gcloud is authenticated."""
    try:
        result = subprocess.run(
            ["gcloud", "auth", "list", "--filter=status:ACTIVE"],
            capture_output=True,
            text=True,
            check=True)
        return "ACTIVE" in result.stdout
    except subprocess.CalledProcessError:
        return False


def validate_gcloud_setup():
    """Validate that gcloud CLI is installed and authenticated.
    
    Raises:
        ValueError: If gcloud CLI is not installed or not authenticated
    """
    if not check_gcloud_installed():
        raise ValueError(
            "gcloud CLI is not installed or not in PATH. "
            "Please install it from: https://cloud.google.com/sdk/docs/install")

    if not check_gcloud_auth():
        raise ValueError("gcloud is not authenticated. "
                         "Please run: gcloud auth login")


def create_gcp_vm(  # pylint: disable=too-many-positional-arguments
    mg_name: str,
    vm_name: str,
    zone: str,
    machine_type: str,
    api_key: str,
    api_url: str,
    spot: bool = True,
    max_idle_time: Optional[Union[int, datetime.timedelta]] = None,
    hostname: Optional[str] = None,
    disk_size: int = 10,
    verbose: bool = True,
):
    """Create a GCP VM instance.
    
    Args:
        mg_name: Name of the machine group
        vm_name: Name of the VM instance
        zone: GCP zone where to create the VM
        machine_type: GCP machine type (e.g., "e2-standard-4")
        api_key: Inductiva API key
        api_url: Inductiva API URL
        spot: Whether to use preemptible (spot) instance
        hostname: Optional hostname for the task runner
        disk_size: Boot disk size in GB, defaults to 10
        verbose: Whether to print verbose output
        max_idle_time: Optional max idle time, defaults to 3 minutes
        
    Raises:
        ValueError: If VM creation fails
    """
    validate_gcloud_setup()

    # Create startup script
    startup_script = create_gcp_startup_script()

    # Write script to temporary file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False) as f:
        f.write(startup_script)
        script_path = f.name

    try:
        if max_idle_time is not None:
            if isinstance(max_idle_time, int):
                max_idle_seconds = max_idle_time * 60
            else:
                max_idle_seconds = int(max_idle_time.total_seconds())
        else:
            max_idle_seconds = 180

        metadata = [
            f"INDUCTIVA_API_KEY={api_key}", f"INDUCTIVA_API_URL={api_url}",
            f"MACHINE_GROUP_NAME={mg_name}",
            f"MAX_IDLE_TIMEOUT={max_idle_seconds}"
        ]

        if hostname:
            metadata.append(f"HOST_NAME={hostname}")

        username = users.get_info().username

        cmd = [
            "gcloud", "compute", "instances", "create", vm_name, "--zone", zone,
            "--machine-type", machine_type, "--image-family", "ubuntu-2204-lts",
            "--image-project", "ubuntu-os-cloud", "--scopes",
            "https://www.googleapis.com/auth/cloud-platform", "--metadata",
            ",".join(metadata), "--metadata-from-file",
            f"startup-script={script_path}", "--boot-disk-size",
            f"{disk_size}GB", "--labels",
            f"managed-by=inductiva,user={username}"
        ]

        if spot:
            cmd.append("--preemptible")

        if verbose:
            print(
                f"Creating GCP VM '{vm_name}' {machine_type} in zone '{zone}'.."
            )

        result = subprocess.run(cmd,
                                capture_output=True,
                                text=True,
                                check=False)

        if result.returncode == 0:
            if verbose:
                print("GCP VM created successfully")
            return
        else:
            raise ValueError(f"Failed to create GCP VM: {result.stderr}")

    except (subprocess.CalledProcessError, OSError) as e:
        raise ValueError(f"Error creating GCP VM: {e}") from e
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass


def delete_gcp_vm(vm_name: str, zone: str, verbose: bool = True):
    """Delete a GCP VM instance.
    
    Args:
        vm_name: Name of the VM instance to delete
        zone: GCP zone where the VM is located
        verbose: Whether to print verbose output
        
    Raises:
        ValueError: If VM deletion fails
    """
    validate_gcloud_setup()

    try:
        cmd = [
            "gcloud", "compute", "instances", "delete", vm_name, "--zone", zone,
            "--quiet"
        ]

        result = subprocess.run(cmd,
                                capture_output=True,
                                text=True,
                                check=False)

        if result.returncode == 0:
            if verbose:
                print(f"GCP VM '{vm_name}' terminated successfully.")
            return
        else:
            # VM might already be terminated
            if "not found" in result.stderr.lower():
                if verbose:
                    print(f"GCP VM '{vm_name}' was already terminated.")
                return
            else:
                raise ValueError(f"Failed to terminate GCP VM: {result.stderr}")

    except (subprocess.CalledProcessError, OSError) as e:
        raise ValueError(f"Error terminating GCP VM: {e}") from e
