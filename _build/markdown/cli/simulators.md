# simulators

## inductiva simulators - CLI interface

```default
inductiva simulators [-h] {list} ...
```

Information about available simulators.

The `inductiva simulators` command provides utility sub-commands for managing the available simulators within the Inductiva API.

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva simulators list (ls)

```default
inductiva simulators list [-h] [--dev]
```

The `inductiva simulators list` sub-command lists all available simulators and associated versions for use with the Inductiva API, including those available for development purposes.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--dev`**]() - include development versions of the simulators.

#### Examples

```bash

$ inductiva simulators list
AVAILABLE SIMULATORS AND VERSIONS FOR PRODUCTION RUNS:

SIMULATOR             VERSIONS
amr-wind              1.4.0, 1.4.0_gpu, 3.4.0, 3.4.0_gpu, 3.4.1, 3.4.1_gpu
apptainer-converter   0.1.0
cans                  2.3.4, 2.4.0, 2.4.0_gpu, 3.0.0, 3.0.0_gpu
cm1                   18, 21.1
coawst                3.8
cp2k                  2025.1, 2025.1_gpu
delft3d               6.04.00
dualsphysics          5.2.1, 5.2.1_gpu, 5.4.1, 5.4.1_gpu
...
```

---
::docsbannersmall
::
