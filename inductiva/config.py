"""Class related to the Web API's connection configuration."""
import inductiva_web_api_client


class Configuration:
    """Class that holds the configuration of the Inductiva WEB API client.

    Attributes:
        api_config: Instance of inductiva_web_api_client's configuration.
        output_dir: Directory to store output files retrieved from the API.
    """

    def __init__(self, address, output_dir):
        self.api_config = inductiva_web_api_client.Configuration(host=address)
        self.output_dir = output_dir
