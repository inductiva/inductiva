**Find answers to commonly asked questions about AMR-Wind.**

<br>

# FAQ

## 1. How can we open `.plt` files generated by AMR-Wind? Is there a Python library available for this?
AMR-Wind, which is based on the AMReX framework, generates output in the form of plotfile directories 
(e.g., `plt00000`, `plt00001`, etc.). While the files use the `.plt` extension, 
they are *not* Tecplot files. These are actually **AMReX plotfiles**, which can be opened 
and processed using the [**yt**](https://yt-project.org/) Python library.

### Example: Loading and working with AMR-Wind plotfiles using `yt`
```python
import yt
import numpy as np

# Load the AMReX plotfile (directory, not a single file)
ds = yt.load("path/to/plt00000")

# Print available fields
print("Available fields:")
for field in ds.field_list:
    print(field)

# Access the full dataset
ad = ds.all_data()

# Extract the 'velocityx' field
velocityx = ad["velocityx"]

# Convert it to a NumPy array matching domain dimensions
velocity_array = np.array(velocityx).reshape(ds.domain_dimensions)

# Extract a 2D xy-slice at z = 2
slice_xy = velocity_array[:, :, 2]
print("XY Slice shape:", slice_xy.shape)
```

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
