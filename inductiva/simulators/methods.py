"""Methods to interact with the available simulators."""
import requests
from collections import defaultdict
import re

# Regex to extract metadata from the image name
# Example:
#   dualsphysics_v5.2.1_dev -> name=dualsphysics, version=5.2.1, isdev=_dev
#   dualsphysics_v5.2.1 -> name=dualsphysics, version=5.2.1, isdev=None
IMAGE_TOKENIZER = re.compile(r"(?P<name>[a-zA-Z-_0-9]+)"
                             r"_v(?P<version>.+?(?=(_dev|$)))"
                             r"(?P<isdev>_dev)?")

DOCKERHUB_URL = "https://hub.docker.com/v2/repositories/inductiva/kutu/tags/"

SKIP_IMAGES = ["base-image", "echo"]


def list_available_images():
    """List available images on DockerHub for the Inductiva API.

    Returns:
        A dictionary with the available images for the API, grouped by branch
        (dev/prod) and simulator name.

    Example:
    >>> list_available_images()
    {
        {'dev':
            {'amr-wind': ['1.4.0'],
              ...
             'dualsphysics': ['5.2.1']
            },
        'prod':
            {'amr-wind': ['1.3.0', '1.4.0'],
              ...
             'dualsphysics': ['5.2.1']
            }
        }
    """

    images = defaultdict(lambda: defaultdict(list))
    url = f"{DOCKERHUB_URL}?page_size=100"

    while url is not None:
        resp = requests.get(url, timeout=10)
        data = resp.json()

        for tag in data["results"]:
            if (match := IMAGE_TOKENIZER.match(tag["name"])) is None:
                continue

            if match.group("name") in SKIP_IMAGES:
                continue

            name, version = match.group("name"), match.group("version")
            branch = "development" if match.group("isdev") else "production"
            images[branch][name].append(version)

        url = data["next"]

    return images
