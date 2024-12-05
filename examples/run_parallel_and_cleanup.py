"""Script to run all example files"""
import os
import subprocess
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

#checks if we are in the examples directory
if not os.path.basename(os.getcwd()) == "examples":
    print("Run this script from the examples directory.")
    sys.exit(1)
# Parse command-line arguments
if len(sys.argv) > 2:
    print("Usage: python run_parallel_and_cleanup.py [max_threads]")
    sys.exit(1)

# Default to 2 threads if no argument is provided
max_threads = int(sys.argv[1]) if len(sys.argv) == 2 else 4
if max_threads <= 0:
    print("Maximum number of threads must be a positive integer.")
    sys.exit(1)

# Define the log directory
log_dir = "logs"
os.makedirs(log_dir, exist_ok=True)


# Function to run a Python script and log its output
def run_script(file_path):
    script_name = os.path.basename(file_path)
    log_file = os.path.join(log_dir, f"{script_name}.log")
    print(f"Running {file_path}...")

    try:
        # Run the script and capture stdout/stderr
        result = subprocess.run(["python3", file_path],
                                capture_output=True,
                                text=True,
                                check=True)
        # Write logs to the log file
        with open(log_file, "w", encoding="utf-8") as log:
            log.write(f"--- Output for {script_name} ---\n")
            log.write("STDOUT:\n")
            log.write(result.stdout + "\n")
            log.write("STDERR:\n")
            log.write(result.stderr + "\n")
        print(f"Logs for {script_name} saved to {log_file}")
    except subprocess.CalledProcessError as e:
        # Handle errors and write to log
        with open(log_file, "w", encoding="utf-8") as log:
            log.write(f"--- Error running {script_name} ---\n")
            log.write("STDOUT:\n")
            log.write(e.stdout + "\n")
            log.write("STDERR:\n")
            log.write(e.stderr + "\n")
        print(f"Error logs for {script_name} saved to {log_file}")


# Step 1: Find all Python files, excluding this script
current_script = os.path.basename(__file__)
python_files = [
    os.path.join(root, file)
    for root, _, files in os.walk(".")
    for file in files
    if file.endswith(".py") and file != current_script
]

# Step 2: Run all scripts in parallel with a limit on active threads
with ThreadPoolExecutor(max_workers=max_threads) as executor:
    futures = [executor.submit(run_script, file) for file in python_files]
    for future in as_completed(futures):
        future.result()  # Wait for each script to complete

# Step 3: Delete folders matching "*input-example"
for root, dirs, _ in os.walk(".", topdown=False):
    for dir_name in dirs:
        if dir_name.endswith("input-example"):
            dir_path = os.path.join(root, dir_name)
            print(f"Deleting directory {dir_path}...")
            shutil.rmtree(dir_path)

# Step 4: Delete inductiva_output folder
shutil.rmtree("inductiva_output")

print("All Python files have been executed with a thread limit of"
      f" {max_threads}, and matching directories have been deleted.")
