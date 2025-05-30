"""Launches a Task-Runner on aws via CLI."""

import argparse
import datetime
import json
import os
import sys
from string import Template
from typing import TextIO

import boto3
import botocore

from inductiva import _api_key, _cli


def _fetch_ubuntu_ami_id(region, profile):
    """Fetches the latest Ubuntu 24.04 LTS AMI ID for a particular region."""
    try:
        session = boto3.Session(profile_name=profile, region_name=region)
        ssm = session.client("ssm")
        parameter = ssm.get_parameter(
            Name=("/aws/service/canonical/ubuntu/server/24.04/stable/"
                  "current/amd64/hvm/ebs-gp3/ami-id"))
        return parameter["Parameter"]["Value"]
    except botocore.exceptions.ClientError as error:
        print(f"Error fetching AMI ID: {error}")
        return None


def _configure_default_security_group(region, profile):
    """Fetches the security group and adds rules for SSH, HTTP, and HTTPS."""
    session = boto3.Session(profile_name=profile, region_name=region)
    ec2_client = session.client("ec2")
    try:
        response = ec2_client.describe_security_groups(Filters=[{
            "Name": "group-name",
            "Values": ["default"]
        }])
        security_group_id = response["SecurityGroups"][0]["GroupId"]
        # Authorize SSH (22), HTTP (80), and HTTPS (443) ports
        ec2_client.authorize_security_group_ingress(
            GroupId=security_group_id,
            IpPermissions=[
                {
                    "IpProtocol": "tcp",
                    "FromPort": 22,
                    "ToPort": 22,
                    "IpRanges": [{
                        "CidrIp": "0.0.0.0/0"
                    }],
                },
                {
                    "IpProtocol": "tcp",
                    "FromPort": 80,
                    "ToPort": 80,
                    "IpRanges": [{
                        "CidrIp": "0.0.0.0/0"
                    }],
                },
                {
                    "IpProtocol": "tcp",
                    "FromPort": 443,
                    "ToPort": 443,
                    "IpRanges": [{
                        "CidrIp": "0.0.0.0/0"
                    }],
                },
            ],
        )
        return security_group_id
    except botocore.exceptions.ClientError as error:
        if "InvalidPermission.Duplicate" in str(error):
            pass
        else:
            print(f"Error configuring security group: {error}")
        return security_group_id


def _create_key_pair(region, profile, key_format):
    """Creates a new RSA key pair in the specified format (PEM or PPK)."""
    session = boto3.Session(profile_name=profile, region_name=region)
    ec2_client = session.client("ec2")
    key_name = "KeyPair-" + datetime.datetime.now().strftime(
        "%Y_%m_%d_%H_%M_%S")
    if key_format not in ["pem", "ppk"]:
        raise ValueError("Invalid key format. Choose 'pem' or 'ppk'.")
    try:
        key_pair = ec2_client.create_key_pair(
            KeyName=key_name,
            KeyType="rsa",
            KeyFormat=key_format,  # PEM or PPK
        )
        # Set the file extension and path based on the format
        file_extension = "pem" if key_format == "pem" else "ppk"
        key_file = os.path.expanduser(f"~/{key_name}.{file_extension}")
        with open(key_file, "w", encoding="UTF8") as file:
            file.write(key_pair["KeyMaterial"])
        os.chmod(key_file, 0o400)
        return key_name
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error creating {key_format.upper()} key pair: {e}")
        return None


def create_iam_role_with_admin_access(role_name, profile):
    session = boto3.Session(profile_name=profile)
    iam_client = session.client("iam")
    try:
        iam_client.get_role(RoleName=role_name)
    except iam_client.exceptions.NoSuchEntityException:
        print(f"Creating role {role_name} with AdministratorAccess...")
        trust_policy = {
            "Version":
                "2012-10-17",
            "Statement": [{
                "Effect": "Allow",
                "Principal": {
                    "Service": "ec2.amazonaws.com"
                },
                "Action": "sts:AssumeRole",
            }],
        }
        iam_client.create_role(
            RoleName=role_name,
            AssumeRolePolicyDocument=json.dumps(trust_policy),
            Description="Role for EC2 instances with Administrator Access",
        )
        iam_client.attach_role_policy(
            RoleName=role_name,
            PolicyArn="arn:aws:iam::aws:policy/AdministratorAccess",
        )
    instance_profile_name = f"{role_name}-instance-profile"
    try:
        iam_client.get_instance_profile(
            InstanceProfileName=instance_profile_name)
    except iam_client.exceptions.NoSuchEntityException:
        iam_client.create_instance_profile(
            InstanceProfileName=instance_profile_name)
        iam_client.add_role_to_instance_profile(
            InstanceProfileName=instance_profile_name, RoleName=role_name)
    return instance_profile_name


def launch_task_runner(args):
    ami_id = _fetch_ubuntu_ami_id(args.region, args.profile)
    security_group_id = _configure_default_security_group(
        args.region, args.profile)
    key_name = _create_key_pair(args.region, args.profile, args.key_format)
    instance_profile_name = create_iam_role_with_admin_access(
        "EC2RoleAdmin", args.profile)

    session = boto3.Session(profile_name=args.profile, region_name=args.region)
    ec2_client = session.client("ec2")

    current_file_path = os.path.dirname(__file__)
    user_data_path = os.path.join(current_file_path, "user_data.sh")

    with open(user_data_path, encoding="UTF8") as f:
        template = Template(f.read())

    user_data_rendered = template.substitute(
        INDUCTIVA_API_KEY=_api_key.get(),
        MACHINE_GROUP_NAME=args.machine_group_name,
    )

    for i in range(args.num_machines):
        unique_name = (f"VM-{i + 1}-"
                       f"{datetime.datetime.now().strftime('%Y_%m_%d_%H_%M')}")
        try:
            response = ec2_client.run_instances(
                ImageId=ami_id,
                InstanceType=args.vm_type,
                KeyName=key_name,
                SecurityGroupIds=[security_group_id],
                UserData=user_data_rendered,
                MinCount=1,
                MaxCount=1,
                IamInstanceProfile={"Name": instance_profile_name},
                TagSpecifications=[{
                    "ResourceType": "instance",
                    "Tags": [{
                        "Key": "Name",
                        "Value": unique_name
                    }],
                }],
                BlockDeviceMappings=[{
                    "DeviceName": "/dev/sda1",
                    "Ebs": {
                        "VolumeSize":
                            args.volume_size,
                        "VolumeType":
                            args.volume_type,
                        "DeleteOnTermination":
                            True,
                        "Iops":
                            args.iops if args.volume_type
                            in ["gp3", "io1", "io2"] else None,
                        "Throughput":
                            args.throughput
                            if args.volume_type == "gp3" else None,
                    },
                }],
            )
            for instance in response["Instances"]:
                print(f"Machine {i + 1} with Name '{unique_name}'"
                      f"and ID '{instance['InstanceId']}' launched at AWS.")
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error launching VM: {e}")


def register(parser):
    """Register the launch task-runner command."""
    subparser = parser.add_parser(
        "launch-aws",
        help="Launches a Task-Runner on aws.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    subparser.description = (
        "The `inductiva task-runner launch-aws` command provides "
        "a way to launch a Task-Runner on aws EC2.")

    _cli.utils.add_watch_argument(subparser)
    subparser.add_argument(
        "machine_group_name",
        type=str,
        help="Name of the machine group to launch the Task-Runner.",
    )

    subparser.add_argument(
        "--vm_type",
        type=str,
        default="t2.micro",
        help="Virtual Machine Type [default: t2.micro].",
    )
    subparser.add_argument(
        "--volume_size",
        type=int,
        default=8,
        help="Size of the root volume in GiB [default: 8]",
    )
    subparser.add_argument(
        "--volume_type",
        type=str,
        choices=["gp3", "gp2", "io1", "io2"],
        default="gp3",
        help="Type of the root volume [default: gp3]",
    )
    subparser.add_argument(
        "--iops",
        type=int,
        default=3000,
        help=
        "IOPS for the volume (applies to io1, io2, gp3 only) [default: 3000]",
    )
    subparser.add_argument(
        "--throughput",
        type=int,
        default=125,
        help="Throughput in MiB/s (gp3 only) [default: 125]",
    )
    subparser.add_argument(
        "--num_machines",
        type=int,
        default=1,
        help="Number of Virtual Machines to launch [default: 1]",
    )
    subparser.add_argument(
        "--region",
        type=str,
        default="eu-west-2",
        help="AWS Region [default: eu-west-2]",
    )
    subparser.add_argument(
        "--key_format",
        type=str,
        choices=["pem", "ppk"],
        default="pem",
        help="Key Pair Format: pem [default] or ppk",
    )
    subparser.add_argument(
        "--profile",
        type=str,
        default="default",
        help="AWS CLI profile [default: default]",
    )

    subparser.set_defaults(func=launch_task_runner)
