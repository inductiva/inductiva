**Find answers to commonly asked questions about OpenSees.**

<br>

# FAQ

## 1. How do I run the EESD OpenSees distribution?
To run the EESD version of OpenSees, simply set the `interface` parameter to "eesd" and specify the appropriate version for your use case.

```python
# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="eesd",
    version="3.0.2")
```

<br>

## 2. How do I run OpenSeesPy?
To run OpenSees scripts written in Python (OpenSeesPy), set the `interface` parameter to "python" and use a version that supports Python.

```python
# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="python",
    version="3.7.1")
```

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
