# Markdown Test File

## pytest-markdown-docs

### Not python example
This first example won't run because it is not a python snippet.

```c
#include <stdio.h>
int main() {
   // printf() displays the string inside quotation
   printf("Hello, World!");
   return 0;
}
```

### No run example
This snippet is python but we define that it should not be tested using the
`notest` tag

```python notest
task_1.download_outputs()
```

### Normal example

This is the first example that will run.
This snippet uses a fixture, configured in the file `conftest.py`.

```python fixture:machine_group_sample
import inductiva

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fds-input-example.zip", unzip=True)

fds = inductiva.simulators.FDS()

task = fds.run(input_dir=input_dir,
               on=machine_group_sample,
               sim_config_filename="mccaffrey.fds",
               post_processing_filename="mccaffrey.ssf",
               n_vcpus=1)
```

### Dependency test
This snippet is dependent on the last snippet.

```python continuation fixture:machine_group_sample
print("fds task", task.get_status())
```

### Fixture test
This snippet is not dependent on any other snippet. However, for it to be tested
we need to provide it with some fixtures (variables or functions that are not
defined inside the snippet):

```python fixture:simulator fixture:input_dir fixture:machine_group_sample

import inductiva

task = simulator.run(input_dir=input_dir,
               on=machine_group_sample,
               sim_config_filename="mccaffrey.fds",
               post_processing_filename="mccaffrey.ssf",
               n_vcpus=1)
```

## Link checker

Internet link [Inductiva](https://inductiva.ai)
Relative link [Home](./Home.md)
