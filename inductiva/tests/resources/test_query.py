"""Unit tests for the inductiva.resources.query."""
from unittest.mock import patch
import pytest

import inductiva


def fake_get_available_machine_types(provider):
    if provider == "ICE":
        return [{
            "machine_type": "g2-standard-4",
            "num_cpus": 4,
            "ram_gb": 16,
            "price": 0.0,
            "provider_id": "ICE",
            "spot": False,
            "region": "porto"
        }, {
            "machine_type": "c2-highmem-8",
            "num_cpus": 8,
            "ram_gb": 64,
            "price": 0.0,
            "provider_id": "ICE",
            "spot": False,
            "region": "porto"
        }]
    elif provider == "GCP":
        return [{
            "machine_type": "n2d-highcpu-128",
            "num_cpus": 128,
            "ram_gb": 128,
            "price": 4.39168,
            "provider_id": "GCP",
            "spot": False,
            "region": "europe-west1"
        }, {
            "machine_type": "n2d-highcpu-16",
            "num_cpus": 16,
            "ram_gb": 16,
            "price": 0.54896,
            "provider_id": "GCP",
            "spot": False,
            "region": "europe-west1"
        }, {
            "machine_type": "n2d-highcpu-2",
            "num_cpus": 2,
            "ram_gb": 2,
            "price": 0.06862,
            "provider_id": "GCP",
            "spot": False,
            "region": "europe-west1"
        }, {
            "machine_type": "n2d-highmem-4",
            "num_cpus": 4,
            "ram_gb": 32,
            "price": 0.25092,
            "provider_id": "GCP",
            "spot": False,
            "region": "europe-west1"
        }, {
            "machine_type": "n2d-highmem-48",
            "num_cpus": 48,
            "ram_gb": 384,
            "price": 3.01104,
            "provider_id": "GCP",
            "spot": False,
            "region": "europe-west1"
        }, {
            "machine_type": "c2d-highcpu-112",
            "num_cpus": 112,
            "ram_gb": 224,
            "price": 0.81424,
            "provider_id": "GCP",
            "spot": True,
            "region": "europe-west1"
        }, {
            "machine_type": "c2d-highcpu-16",
            "num_cpus": 16,
            "ram_gb": 32,
            "price": 0.11632,
            "provider_id": "GCP",
            "spot": True,
            "region": "europe-west1"
        }, {
            "machine_type": "c2d-highcpu-2",
            "num_cpus": 2,
            "ram_gb": 4,
            "price": 0.01454,
            "provider_id": "GCP",
            "spot": True,
            "region": "europe-west1"
        }]


@patch("inductiva.resources.methods.get_available_machine_types",
       new=fake_get_available_machine_types)
def test_resources__query__ice():
    res = inductiva.resources.query(provider="ICE")
    assert len(res) == 2
    assert res[0]["machine_type"] == "g2-standard-4"
    assert res[1]["machine_type"] == "c2-highmem-8"


@patch("inductiva.resources.methods.get_available_machine_types",
       new=fake_get_available_machine_types)
def test_resources__query__gcp():
    res = inductiva.resources.query(provider="GCP")
    assert len(res) == 8
    assert res[0]["machine_type"] == "n2d-highcpu-128"
    assert res[1]["machine_type"] == "n2d-highcpu-16"
    assert res[2]["machine_type"] == "n2d-highcpu-2"


@patch("inductiva.resources.methods.get_available_machine_types",
       new=fake_get_available_machine_types)
def test_resources__query__gcp_string_filter():
    with pytest.raises(ValueError):
        inductiva.resources.query(provider="GCP", query_filter="all")


def c2d_filter(machine_type):
    return "c2d" in machine_type


@patch("inductiva.resources.methods.get_available_machine_types",
       new=fake_get_available_machine_types)
def test_resources__query__gcp_c2d():
    res = inductiva.resources.query(provider="GCP", query_filter=c2d_filter)
    assert len(res) == 3
    assert res[0]["machine_type"] == "c2d-highcpu-112"
    assert res[1]["machine_type"] == "c2d-highcpu-16"
    assert res[2]["machine_type"] == "c2d-highcpu-2"
