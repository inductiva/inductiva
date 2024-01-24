# API access tokens

To use Inductiva API you will need a valid API token. At the moment, you can request an API token [here](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform). 

Once you have your API token, you can either:
- set the `INDUCTIVA_API_KEY` as an environment variable in your terminal as follows

```bash
export INDUCTIVA_API_KEY="YOUR_API_KEY"
```

- set the `api_key` attribute of the `inductiva` module in your Python script as follows

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"
```

In the latter case, you will need to set your API key in all the Python scripts that use Inductiva API.
Therefore, we recommend the first approach and that is the one we assume in the code snippets of this documentation.

You are good to go! You can start exploring Inductiva API with the examples below.