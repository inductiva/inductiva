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
    script_content = """#!/bin/bash

set -e

# -------------------------------
# 1. Install Docker
# -------------------------------
apt-get update -y
apt-get install -y docker.io

systemctl enable docker
systemctl start docker

# -------------------------------
# 2. Prepare directories and volume
# -------------------------------
mkdir -p /home/runner/apptainer
chmod 777 /home/runner/apptainer

docker volume create workdir

# -------------------------------
# 3. Pull Docker images
# -------------------------------
docker pull inductiva/task-runner:latest
docker pull inductiva/file-tracker:latest

# -------------------------------
# 4. Read metadata and export env variables
# -------------------------------
for var in INDUCTIVA_API_KEY INDUCTIVA_API_URL MACHINE_GROUP_NAME; do
    export "$var"=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$var" -H "Metadata-Flavor: Google")
    echo "$var=${!var}"
done

# -------------------------------
# 5. Launch file-tracker
# -------------------------------
docker run -d --name file-tracker \\
  --network host \\
  -v workdir:/workdir \\
  -e API_URL="$INDUCTIVA_API_URL" \\
  -e USER_API_KEY="$INDUCTIVA_API_KEY" \\
  inductiva/file-tracker:latest

# -------------------------------
# 6. Launch task-runner
# -------------------------------
docker run -d --name task-runner \\
  --network host \\
  --privileged \\
  --platform linux/amd64 \\
  -v /home/runner/apptainer:/executer-images \\
  -v workdir:/workdir \\
  --add-host host.docker.internal:host-gateway \\
  -e EXECUTER_IMAGES_DIR=/executer-images \\
  -e API_URL="$INDUCTIVA_API_URL" \\
  -e USER_API_KEY="$INDUCTIVA_API_KEY" \\
  -e MACHINE_GROUP_NAME="$MACHINE_GROUP_NAME" \\
  -e HOST_NAME="${TASK_RUNNER_HOSTNAME:-$(hostname)}" \\
  inductiva/task-runner:latest

# -------------------------------
# 7. Monitor task-runner and shutdown VM if it stops
# -------------------------------
METADATA_URL="http://metadata.google.internal/computeMetadata/v1/instance"
VM_NAME=$(curl -s -H "Metadata-Flavor: Google" "$METADATA_URL/name")
ZONE=$(curl -s -H "Metadata-Flavor: Google" "$METADATA_URL/zone" | awk -F/ '{print $4}')

while true; do
    sleep 10
    if [ -z "$(docker ps -q -f name=task-runner)" ]; then
        echo "Task runner stopped, deleting VM..."
        gcloud compute instances delete "$VM_NAME" --zone="$ZONE" --quiet
        break
    fi
done
"""
    return script_content


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
            metadata.append(f"TASK_RUNNER_HOSTNAME={hostname}")

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
            print(f"Creating GCP VM '{vm_name}' in zone '{zone}'...")
            print(f"Machine type: {machine_type}")
            if spot:
                print("Using spot instance")

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
