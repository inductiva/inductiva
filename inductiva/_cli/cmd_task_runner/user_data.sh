#!/bin/bash

# Enable unprivileged user namespaces (required for Apptainer/Singularity containers)
sysctl -w kernel.unprivileged_userns_clone=1
echo 'kernel.unprivileged_userns_clone=1' > /etc/sysctl.d/99-unpriv-userns.conf
sysctl --system

# Update system and install packages
apt-get update -y
apt-get install -y ca-certificates curl python-is-python3 python3-pip python3-venv

# Install Docker
install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
chmod a+r /etc/apt/keyrings/docker.asc

echo "deb [arch=$$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] \
https://download.docker.com/linux/ubuntu \
$$(. /etc/os-release && echo "$$VERSION_CODENAME") stable" > /etc/apt/sources.list.d/docker.list

apt-get update -y
apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

usermod -aG docker ubuntu

python3 -m venv /opt/inductiva-env
/opt/inductiva-env/bin/pip install --upgrade pip
/opt/inductiva-env/bin/pip install inductiva[task-runner]

INDUCTIVA_API_KEY="$INDUCTIVA_API_KEY" /opt/inductiva-env/bin/inductiva task-runner launch "$MACHINE_GROUP_NAME"
