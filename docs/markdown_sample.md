# Markdown Test File

### Not python example
This first example wont run because its not a python snippet.

```c
#include <stdio.h>
int main() {
   // printf() displays the string inside quotation
   printf("Hello, World!");
   return 0;
}
```

### No run example
This snippet is python but we define it not to run using 

```python notest
task_1.download_outputs()
```

### Normal example
This is the first example that will run. This snippet is just the normal case where there are no dependencies or ficture vars.

```python
import inductiva

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fds-input-example.zip", unzip=True)

fds = inductiva.simulators.FDS()

task = fds.run(input_dir=input_dir,
               sim_config_filename="mccaffrey.fds",
               post_processing_filename="mccaffrey.ssf",
               n_vcpus=1)
```

### Dependency test
This snipper is dependant on the last snippet.

```python continuation
print("fds task",task.get_status())
```

### Fixture test
This snippet is not dependant on any other snippet. But, we need to provide it with some fixtures (variables or functions that are not defined inside the snippet).

```python fixture:simulator fixture:input_dir

import inductiva

task = simulator.run(input_dir=input_dir,
               sim_config_filename="mccaffrey.fds",
               post_processing_filename="mccaffrey.ssf",
               n_vcpus=1)
```