"""Script for demonstrating killing a blocking API call.

The script runs a demo task in the API that sleeps for a specified amount of
time. Press ctrl+c while it is running to kill the task on the API.

Usage: python kill.py secs={SLEEP_SECS}
"""
import inductiva

from absl import flags
from absl import app

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://localhost:8080",
                    "Base URL of the Inductiva API.")

flags.DEFINE_integer("secs", 5, "Num. secs to sleep.")


def main(_):
    inductiva.api_url = FLAGS.api_url

    inductiva.test.sleep(secs=FLAGS.secs)


if __name__ == "__main__":
    app.run(main)
