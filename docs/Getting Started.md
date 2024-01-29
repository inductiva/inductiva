# Getting Started
Welcome to Inductiva, a Python API package that simplifies running large-scale 
physical system simulations in the cloud. Just few lines of code unleash its potential, freeing 
you from local hardware limits and the technical hassles of scaling simulations. 

### What We'll Cover:

* [Inductiva API Installation]()
* [What to read next]()
* [Troubleshooting]()
* [Frequently Asked Questions (FAQs)]()
* [Useful Terminology]()

## Inductiva API Installation

To install the latest Inductiva package release on PyPI with pip, just open your 
terminal and run:

```
pip install --upgrade inductiva
```
---

> **Tip:** To ensure smooth installation and operation of our API, we recommend 
keeping your Python package manager, Pip, up-to-date. Visit [Pip's official installation guide](https://pip.pypa.io/en/stable/installation/) for detailed instructions on updating Pip.


Encountering issues? Don’t worry. Head over to [Troubleshooting]() to work it out.

Otherwise, you're all set! You can now decide [what to read next]().

## What to read next

You can start experimenting with Inductiva by choosing from a variety of ready-to-use, 
[open-source simulators]() that are immediately available for use on our default 
machines.

You can set up custom, dedicated machine groups for more complex or [larger-scale simulation]() 
tasks, tapping into the more extensive [computational resources]() of Inductiva.

You can get inspired by how other researchers have used Inductiva API to configure, 
execute, and scale their simulations with minimal coding effort by visiting our 
[example gallery]().

You can troubleshoot installation problems that you might encounter with Inductiva 
API by checking [troubleshooting](#troubleshooting), [FAQs](), or [getting in touch]() 
with our support team directly.

## Troubleshooting

If installing the package fails, you can retry it on a new Python virtual environment. 
A [virtual environment](https://docs.python.org/3/library/venv.html) allows you to 
have a fresh Python environment with isolated dependencies. In your shell, run:

```
python -m venv <venv>
```

In that command, you should replace `<venv>` with the path (*e.g.*, `.venv`) in 
which you would like to create the environment. Then, to activate the environment 
(again, correctly replacing `<venv>`), run:

For `bash`/`zsh`:

```
source <venv>/bin/activate
```

For `cmd.exe` (Windows):

```
<venv>\Scripts\activate.bat
```

For `PowerShell` (Windows):
```
<venv>\Scripts\Activate.ps1
```

After activating the virtual environment, you can install the package as described 
below:

```
pip install --upgrade inductiva
```

## Frequently Asked Questions
This section could be removed or updated later.
## Useful Terminology
This documentation uses the following terms:

https://packaging.python.org/en/latest/glossary/ 
