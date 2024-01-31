# Troubleshooting Guide
In this troubleshooting section, we aim to guide you through resolving issues related to incorrect or incomplete [installation of the Inductiva package]() or its dependencies.
### What We'll Cover:

* [Installation Failures]()
* [Frequently Asked Questions (FAQs)]()
* [Useful Terminology]()

## Installation Failures

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

## Frequently Asked Questions (FAQs)

**What should I do if Inductiva does not support the simulation package I need?**

*If the simulation package you wish to use is not supported by Inductiva, please reach out to our support team for assistance at support@inductiva.ai.*

**How do I resolve authentication problems with the Inductiva API?**

*Authentication issues typically arise from invalid credentials or expired tokens. Double-check your credentials and refresh your API tokens. For further assistance, our support team is ready to help.*

**What can I do if I reach my resource allocation limits?**

*If you've hit your quota for computational resources, consider reviewing your current usage. To request an increase in your resource quotas, please reach out to support@inductiva.ai.*

**The simulation package I want to use isn't supported. What are my options?**

**How can I fix configuration errors in my simulation setup?**

**What should I do about dependency conflicts?**

**How can I optimize the performance of my simulations?**
Performance optimization involves adjusting simulation parameters. Our documentation includes tips for efficiently utilizing [Inductiva's computational resources]().


## Useful Terminology
This documentation uses the following terms:

https://packaging.python.org/en/latest/glossary/ 

