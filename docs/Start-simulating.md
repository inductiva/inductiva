# Start simulating

After having installed and prepared your environment, you can start simulating with Inductiva API.

## Simulate with Inductiva API

Follow the example with:

```python
import inductiva

# Download the configuration files
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True
)

# Initialize the Simulator
simulator = inductiva.simulators.REEF3D()

# Run the simulation
task = simulator.run(input_dir=input_dir)
task.wait()
```
