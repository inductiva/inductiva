"""GCP BYOC (Bring Your Own Cloud) utilities."""
import subprocess
import tempfile
import os
import datetime
import time
import logging
import uuid
import re

import inductiva
from inductiva.utils import format_utils


def register(machine_group):
    """Register machine group for client-side GCP management."""
    logging.info("■ Registering %s configurations (client-side):",
                 machine_group.short_name())

    machine_group._id = str(uuid.uuid4())
    machine_group._name = f"client-mg-{machine_group.machine_type}-{machine_group._id[:8]}"

    _client_vm_info = {
        'vm_name': machine_group._name,
        'zone': machine_group.zone,
        'status': 'registered',
        'created_at': datetime.datetime.now()
    }
    machine_group._client_vm_info = _client_vm_info

    machine_group.create_time = datetime.datetime.now()
    machine_group.num_machines = getattr(machine_group, 'num_machines', 1)
    machine_group._active_machines = 0

    if isinstance(machine_group.max_idle_time, int):
        machine_group.max_idle_time = datetime.timedelta(
            minutes=machine_group.max_idle_time)

    machine_group._cpu_info = type(
        'CPUInfo', (), {
            'cpu_cores_logical':
                _estimate_vcpus_from_machine_type(machine_group.machine_type),
            'cpu_cores_physical':
                _estimate_vcpus_from_machine_type(machine_group.machine_type) //
                2
        })()

    machine_group._gpu_info = type('GPUInfo', (), {
        'gpu_count': 0,
        'gpu_name': None
    })()

    machine_group._cost_per_hour = type(
        'CostInfo', (), {
            'min': 0.1,
            'max': 0.1,
            'min_reason': 'Client-side GCP',
            'max_reason': 'Client-side GCP'
        })()

    _register_machine_group_backend(machine_group)

    machine_group._log_machine_group_info()


def _register_machine_group_backend(machine_group):
    """Register machine group with backend for simulator compatibility."""
    logging.info("■ Registering with backend for simulator compatibility...")

    instance_group_config = inductiva.client.models.RegisterVMGroupRequest(
        machine_type=machine_group.machine_type,
        provider_id="LOCAL",  # Register as LOCAL for backend compatibility
        threads_per_core=machine_group.threads_per_core,
        disk_size_gb=machine_group.data_disk_gb,
        max_idle_time=machine_group._timedelta_to_seconds(
            machine_group.max_idle_time),
        auto_terminate_ts=machine_group._convert_auto_terminate_ts(
            machine_group.auto_terminate_ts),
        dynamic_disk_resize_config=machine_group._dynamic_disk_resize_config(),
        custom_vm_image=machine_group._custom_vm_image,
        zone=machine_group.zone,
        disk_auto_delete=machine_group.auto_delete_disk,
        num_vms=machine_group.num_machines,
        spot=machine_group.spot,
        is_elastic=False,
    )

    body = machine_group._api.register_vm_group(
        register_vm_group_request=instance_group_config,)

    machine_group._update_attributes_from_response(body)

    machine_group._cpu_info = type(
        'CPUInfo', (), {
            'cpu_cores_logical':
                _estimate_vcpus_from_machine_type(machine_group.machine_type),
            'cpu_cores_physical':
                _estimate_vcpus_from_machine_type(machine_group.machine_type) //
                2
        })()

    machine_group._gpu_info = type('GPUInfo', (), {
        'gpu_count': 0,
        'gpu_name': None
    })()

    machine_group._cost_per_hour = type(
        'CostInfo', (), {
            'min': 0.1,
            'max': 0.1,
            'min_reason': 'Client-side GCP',
            'max_reason': 'Client-side GCP'
        })()


def _estimate_vcpus_from_machine_type(machine_type):
    """Estimate vCPUs from machine type string."""
    # Extract number from machine type (e.g., "c2d-standard-8" -> 8)
    match = re.search(r'-(\d+)$', machine_type)
    if match:
        return int(match.group(1))
    raise ValueError(
        f"Could not extract vCPU count from machine type: {machine_type}")


def create_gcp_startup_script() -> str:
    """Create the startup script for GCP VM."""
    script_content = f"""#!/bin/bash

set -e

# -------------------------------
# 1. Install Docker
# -------------------------------
apt update
apt install -y docker.io

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
    echo "$var=${{!var}}"
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
  -e HOST_NAME="${{TASK_RUNNER_HOSTNAME:-$(hostname)}}" \\
  inductiva/task-runner:latest

# -------------------------------
# 7. Monitor task-runner and shutdown VM if it stops
# -------------------------------
METADATA_URL="http://metadata.google.internal/computeMetadata/v1/instance"
VM_NAME=$(curl -s -H "Metadata-Flavor: Google" "$METADATA_URL/name")
ZONE=$(curl -s -H "Metadata-Flavor: Google" "$METADATA_URL/zone" | awk -F/ '{{print $4}}')

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
        subprocess.run(['gcloud', '--version'],
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
            ['gcloud', 'auth', 'list', '--filter=status:ACTIVE'],
            capture_output=True,
            text=True,
            check=True)
        return 'ACTIVE' in result.stdout
    except subprocess.CalledProcessError:
        return False


def start(machine_group, verbose: bool = True):
    """Start GCP VMs using client-side management."""
    if not check_gcloud_installed():
        print("Error: gcloud CLI is not installed or not in PATH.")
        print(
            "Please install gcloud CLI: https://cloud.google.com/sdk/docs/install"
        )
        print("Or install with: pip install 'inductiva[gcp]'")
        return False

    if not check_gcloud_auth():
        print("Error: gcloud is not authenticated.")
        print("Please run 'gcloud auth login' to authenticate.")
        return False

    from inductiva import _api_key, api_url
    api_key = _api_key.get()
    if not api_key:
        print("Error: No API key found. Please set your API key first.")
        return False

    logging.info("Starting %s (client-side GCP)...", repr(machine_group))
    start_time = time.time()

    # Create startup script
    startup_script = create_gcp_startup_script()

    # Write script to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
        f.write(startup_script)
        script_path = f.name

    try:
        cmd = [
            'gcloud', 'compute', 'instances', 'create', machine_group._name,
            '--zone', machine_group.zone, '--machine-type',
            machine_group.machine_type, '--image-family', 'ubuntu-2204-lts',
            '--image-project', 'ubuntu-os-cloud', '--scopes',
            'https://www.googleapis.com/auth/cloud-platform', '--metadata',
            f'INDUCTIVA_API_KEY={api_key},INDUCTIVA_API_URL={api_url},MACHINE_GROUP_NAME={machine_group._name}',
            '--metadata-from-file', f'startup-script={script_path}'
        ]

        if machine_group.spot:
            cmd.append('--preemptible')

        if verbose:
            print(
                f"Creating GCP VM '{machine_group._name}' in zone '{machine_group.zone}'..."
            )
            print(f"Machine type: {machine_group.machine_type}")
            if machine_group.spot:
                print("Using preemptible instance (spot pricing)")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            machine_group._started = True
            machine_group._active_machines = machine_group.num_machines
            machine_group._client_vm_info['status'] = 'running'
            machine_group._client_vm_info['started_at'] = datetime.datetime.now(
            )

            creation_time = format_utils.seconds_formatter(time.time() -
                                                           start_time)
            logging.info("%s successfully started in %s.", machine_group,
                         creation_time)

            if verbose:
                print("GCP VM created successfully!")
                print(f"VM Name: {machine_group._name}")
                print(f"Zone: {machine_group.zone}")
                print(
                    "The task-runner will start automatically once the VM is ready."
                )

            return True
        else:
            print("Failed to create GCP VM:")
            print(result.stderr)
            return False

    except Exception as e:
        print(f"Error creating GCP VM: {e}")
        return False
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass


def terminate(machine_group, verbose: bool = True):
    """Terminate GCP VMs using client-side management."""
    if not check_gcloud_installed():
        print("Error: gcloud CLI is not installed or not in PATH.")
        return False

    logging.info("Terminating %s (client-side GCP)...", repr(machine_group))

    try:
        cmd = [
            'gcloud', 'compute', 'instances', 'delete', machine_group._name,
            '--zone', machine_group.zone, '--quiet'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            # Also notify backend that machine group is terminated
            try:
                machine_group._api.delete_vm_group(
                    machine_group_id=machine_group.id)
            except Exception as e:
                logging.warning(f"Failed to notify backend of termination: {e}")

            machine_group._started = False
            machine_group._active_machines = 0
            machine_group._client_vm_info['status'] = 'terminated'
            machine_group._client_vm_info[
                'terminated_at'] = datetime.datetime.now()

            logging.info("%s terminated.", machine_group)

            if verbose:
                print(
                    f"GCP VM '{machine_group._name}' terminated successfully.")

            return True
        else:
            # VM might already be terminated (e.g., by startup script)
            if "not found" in result.stderr.lower():
                # Also notify backend that machine group is terminated
                try:
                    machine_group._api.delete_vm_group(
                        machine_group_id=machine_group.id)
                except Exception as e:
                    logging.warning(
                        f"Failed to notify backend of termination: {e}")

                machine_group._started = False
                machine_group._active_machines = 0
                machine_group._client_vm_info['status'] = 'terminated'
                logging.info("%s was already terminated.", machine_group)
                if verbose:
                    print(
                        f"GCP VM '{machine_group._name}' was already terminated."
                    )
                return True
            else:
                print("Failed to terminate GCP VM:")
                print(result.stderr)
                return False

    except Exception as e:
        print(f"Error terminating GCP VM: {e}")
        return False
