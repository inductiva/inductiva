"""GCP BYOC (Bring Your Own Cloud) utilities."""
import subprocess
import re
import tempfile
import os
from typing import Optional, Tuple


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


def create_gcp_vm(  # pylint: disable=too-many-positional-arguments
        vm_name: str,
        zone: str,
        machine_type: str,
        api_key: str,
        api_url: str,
        spot: bool = True,
        hostname: Optional[str] = None,
        verbose: bool = True) -> Tuple[bool, Optional[str]]:
    """Create a GCP VM instance.
    
    Args:
        vm_name: Name of the VM instance
        zone: GCP zone where to create the VM
        machine_type: GCP machine type (e.g., "e2-standard-4")
        api_key: Inductiva API key
        api_url: Inductiva API URL
        spot: Whether to use preemptible (spot) instance
        hostname: Optional hostname for the task runner
        verbose: Whether to print verbose output
        
    Returns:
        Tuple of (success: bool, error_message: Optional[str])
    """
    if not check_gcloud_installed():
        return False, "gcloud CLI is not installed or not in PATH"

    if not check_gcloud_auth():
        return False, "gcloud is not authenticated"

    # Create startup script
    startup_script = create_gcp_startup_script()

    # Write script to temporary file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False) as f:
        f.write(startup_script)
        script_path = f.name

    try:
        metadata = [
            f"INDUCTIVA_API_KEY={api_key}", f"INDUCTIVA_API_URL={api_url}",
            f"MACHINE_GROUP_NAME={vm_name}"
        ]

        if hostname:
            metadata.append(f"HOST_NAME={hostname}")

        cmd = [
            "gcloud", "compute", "instances", "create", vm_name, "--zone", zone,
            "--machine-type", machine_type, "--image-family", "ubuntu-2204-lts",
            "--image-project", "ubuntu-os-cloud", "--scopes",
            "https://www.googleapis.com/auth/cloud-platform", "--metadata",
            ",".join(metadata), "--metadata-from-file",
            f"startup-script={script_path}"
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
                print("GCP VM created successfully!")
                print(f"VM Name: {vm_name}")
                print(f"Zone: {zone}")
                print("The task-runner will start automatically once the VM "
                      "is ready.")
            return True, None
        else:
            return False, result.stderr

    except (subprocess.CalledProcessError, OSError) as e:
        return False, f"Error creating GCP VM: {e}"
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass


def delete_gcp_vm(vm_name: str,
                  zone: str,
                  verbose: bool = True) -> Tuple[bool, Optional[str]]:
    """Delete a GCP VM instance.
    
    Args:
        vm_name: Name of the VM instance to delete
        zone: GCP zone where the VM is located
        verbose: Whether to print verbose output
        
    Returns:
        Tuple of (success: bool, error_message: Optional[str])
    """
    if not check_gcloud_installed():
        return False, "gcloud CLI is not installed or not in PATH"

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
            return True, None
        else:
            # VM might already be terminated
            if "not found" in result.stderr.lower():
                if verbose:
                    print(f"GCP VM '{vm_name}' was already terminated.")
                return True, None
            else:
                return False, result.stderr

    except (subprocess.CalledProcessError, OSError) as e:
        return False, f"Error terminating GCP VM: {e}"
