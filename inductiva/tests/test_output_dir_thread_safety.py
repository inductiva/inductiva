"""Tests thread safety"""
import os
import time
from concurrent.futures import ThreadPoolExecutor

import inductiva

CURRENT_DIR = os.getcwd()


def thread_1():
    inductiva.set_output_dir("thread_1")
    time.sleep(5)
    inductiva.set_output_dir("thread_1")
    return inductiva.get_output_dir()


def thread_2():
    time.sleep(1)
    inductiva.set_output_dir("thread_2")
    return inductiva.get_output_dir()


def test_thread_safety():
    """Tests that threads do not make definite changes to context."""
    inductiva.set_output_dir("test_dir")
    output_dir = inductiva.get_output_dir()

    with ThreadPoolExecutor(max_workers=2) as executer:
        future_1 = executer.submit(thread_1)
        future_2 = executer.submit(thread_2)

        result_1 = future_1.result()
        result_2 = future_2.result()
    assert inductiva.get_output_dir() == output_dir
    assert result_1 == "thread_1"
    assert result_2 == "thread_2"

    os.chdir(CURRENT_DIR)
