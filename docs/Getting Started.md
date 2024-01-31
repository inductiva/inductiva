# Getting Started with Inductiva
Welcome to Inductiva, a Python API package that enables running popular [open-source physical system simulation packages]() 
in the Cloud. 

This tutorial provides a step-by-step guide to help you get started with Inductiva 
before running your first simulation project. Should you encounter any errors or 
issues along the way, our [troubleshooting guide]() is readily available to assist 
you in resolving them efficiently.

### Steps We'll Cover:

1. [Install the Inductiva API]()
2. [Request Your API Access Token]()
3. [Verify the Installation with a Test Run]()
4. [Examine Project Example]()
5. [Explore what to read next]()

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

## Project Example
Finally, check out this Project Example to demonstrate how you can run simulations 
through the Inductiva API.

## What to read next

Learn how to [run one of your own simulation projects]() for the first time through 
the Inductiva API, on the simulator of your choice.

Learn how to [customize the hardware setup]() you use for running your simulations 
through the Inductiva API, and explore the available hardware options to enhance your project's performance.

If you're looking for inspiration, learn how a group of coastal engineering researchers 
at the University of Porto's Faculty of Engineering have [used the Inductiva API to simulate the most optimal breakwater](https://inductiva.ai/blog/article/scaling-coastal-engineering-projects-inductiva-api) 
to protect some of Portugal's most endangered coastlines.

Troubleshoot installation problems that you might encounter with Inductiva 
API by checking our [troubleshooting guide](#troubleshootingguide), [FAQs](), or [getting in touch]() 
with our support team directly.

