"""
Class related to the Web API's connection configuration.
"""
import inductiva_web_api_client


class Configuration:
    """Class that holds the configuration of the Inductiva WEB API client."""

    def __init__(self, address, output_dir):
        self.api_config = inductiva_web_api_client.Configuration(host=address)
        self.output_dir = output_dir
