"""Script to run all example files"""
import os
import sys
import shutil
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_script(file_path, log_dir="logs"):
    """Run a Python script and log its output to a file."""
    script_name = os.path.basename(file_path)
    log_file = os.path.join(log_dir, f"{script_name}.log")
    print(f"Running {file_path}...")

    stdout = ""
    stderr = ""
    success = False

    try:
        # Run the script and capture stdout/stderr
        result = subprocess.run(["python3", file_path],
                                capture_output=True,
                                text=True,
                                check=True)
        stdout = result.stdout
        stderr = result.stderr
        success = True
    except subprocess.CalledProcessError as e:
        stdout = e.stdout
        stderr = e.stderr
        success = False
    finally:
        # Write logs
        with open(log_file, "w", encoding="utf-8") as log:
            log.write(f"--- Output for {script_name} ---\n")
            log.write("STDOUT:\n")
            log.write(stdout + "\n")
            log.write("STDERR:\n")
            log.write(stderr + "\n")

        status = "\033[92m✔" if success else "\033[91m✘"

        print(f"{status} Logs for {script_name} saved to {log_file}\033[0m")


def gather_python_files(path):
    """Find all Python files in the given directory and its subdirectories.
    
    Excludes the current script from the list of files.
    """
    current_script = os.path.basename(__file__)
    python_files = [
        os.path.join(root, file)
        for root, _, files in os.walk(path)
        for file in files
        if file.endswith(".py") and file != current_script
    ]
    return python_files


def run_python_files(python_files, max_working_threads):
    """Run all Python files in parallel with a limit on active threads."""
    with ThreadPoolExecutor(max_workers=max_working_threads) as executor:
        _ = list(executor.map(run_script, python_files))


def delete_input_examples_from_folder(path):
    """Delete all directories named 'input-example' from the given path."""
    for root, dirs, _ in os.walk(path, topdown=False):
        for dir_name in dirs:
            if dir_name.endswith("input-example"):
                dir_path = os.path.join(root, dir_name)
                print(f"Deleting directory {dir_path}...")
                shutil.rmtree(dir_path)


def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run scripts in parallel and clean up.")
    parser.add_argument("max_threads",
                        nargs="?",
                        type=int,
                        help="Maximum number of threads to use (default: 4)")

    args = parser.parse_args()

    #checks if we are in the examples directory
    if not os.path.basename(os.getcwd()) == "examples":
        print("Run this script from the examples directory.")
        sys.exit(1)

    max_threads = args.max_threads
    if max_threads is not None and max_threads <= 0:
        print("Maximum number of threads must be a positive integer.")
        sys.exit(1)

    # Define the log directory
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)

    all_python_files = gather_python_files(".")

    run_python_files(all_python_files, max_threads)

    delete_input_examples_from_folder(".")

    shutil.rmtree("inductiva_output")

    print("All Python files have been executed with a thread limit of"
          f" {max_threads}, and matching directories have been deleted.")

    #write this in red
    print("\033[91m"
          "IMPORTANT: This script only runs the examples.\n"
          "It does not check if the examples ended with success or failure.\n"
          "Please check your tasks manually to ensure they completed "
          "successfully.\033[0m")


if __name__ == "__main__":
    sys.exit(main())
