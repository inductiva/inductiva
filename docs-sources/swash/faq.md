**Find answers to commonly asked questions about SWAN.**

<br>

# FAQ

## 1. My simulation runs almost to the end, but then fails and the status changes to “Failed.” Why?
One of the most common reasons for this issue is a mismatch between the version of SWASH you are requesting Inductiva to run (e.g., `11.01`) and the version for which your configuration file was originally created (e.g., `9.01`).

Make sure you explicitly specify the correct version when creating the SWASH object. For example:

```python
# Initialize the Simulator
swash = inductiva.simulators.SWASH(version="10.05")
```

You can check the [list of currently available versions](versions-and-containers.md) to ensure compatibility.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
