import os
from concurrent.futures import ThreadPoolExecutor

import inductiva


def test_thread_safety():
    """Tests that threads do not make definite changes to context."""
    inductiva.set_output_dir("test_dir")
    output_dir = inductiva.get_output_dir()

    num_threads = 5
    output_dirs = [f"output_dir_{i}" for i in range(num_threads)]

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [
            executor.submit(inductiva.set_output_dir, output_dirs[i])
            for i in range(num_threads)
        ]

        for future in futures:
            future.result()

    # Test that the output_dir remained unchanged.
    assert inductiva.get_output_dir() == output_dir
