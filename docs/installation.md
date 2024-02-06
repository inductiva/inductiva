# Installation

Welcome to Inductiva, a Python API package that enables running popular
[open-source physical system simulation packages]() in the Cloud. 

This is a step-by-step guide to help you get started with Inductiva 
before running your first simulation project. Should you encounter any errors or 
issues along the way, our [troubleshooting guide]() is readily available to assist 
you in resolving them efficiently.

## Install the Inductiva API with pip
If you've set up the API access token, you're ready to install the 
latest Inductiva package release on PyPI with pip. 

Simply open your terminal and run:

```
pip install --upgrade inductiva
```

Encountering issues? Don’t worry. Head over to our [troubleshooting guide]() to 
work it out.

## Request Your API Access Token
If you don't have a valid API access token, [fill this form](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform) and request 
your own personal API key.

Once you receive your API token, you can then set the `INDUCTIVA_API_KEY` as an 
environment variable in your terminal:
```
export INDUCTIVA_API_KEY="YOUR_API_KEY"
```

## Verify the Installation with a Test Run

Next, make sure everything is set up correctly with a quick test run.
You have two options to test the installation:
1. through the command line interface (CLI) tool that gets installed
   automatically when the `pip` command runs;
2. programmatically via direct import of the `inductiva` package in a python
   command.

Both options are shown below and both output the same version:

```console
# using the CLI tool ($ is the shell prompt):
$ inductiva --version
inductiva 0.4.2
# programmatically through direct python invocation:
$ python -c 'import inductiva; print(inductiva.__version__)'
0.4.2
```


Troubleshoot installation problems that you might encounter with Inductiva 
API by checking our [troubleshooting guide](#troubleshootingguide), [FAQs](), or
[getting in touch]() with our support team directly.

