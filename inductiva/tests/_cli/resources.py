"""Tests for cli resources commands"""
from unittest import mock
import pytest
import subprocess

import inductiva

def test():
    res = subprocess.run(["inductiva", "resources","terminate","abc"], stdout=subprocess.PIPE)
    
    assert "No active" in res.stdout.decode("utf-8")