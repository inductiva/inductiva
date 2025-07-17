# inductiva **simulators** [\[subcommands\]](#subcommands) [\[flags\]](#flags)

The `inductiva simulators` command provides utility subcommands 
for managing the available simulators within the Inductiva API.

````{eval-rst}
.. seealso::
   For complete API documentation, see the `Simulators <https://inductiva.ai/guides/api-functions/api/inductiva.simulators>`_ class documentation
````

## Subcommands
### `list (ls)`
List the available simulators in the Inductiva API, including the supported versions.

```bash
inductiva simulators list
```

or using the alias:

```bash
inductiva simulators ls
```

Sample output:

```bash
$ inductiva simulators list
AVAILABLE SIMULATORS AND VERSIONS FOR PRODUCTION RUNS:

 SIMULATOR             VERSIONS
 amr-wind              3.4.0, 1.4.0
 apptainer-converter   0.1.0
 cans                  2.4.0, 2.3.4
 cm1                   21.1, 18
 coawst                3.8
 cp2k                  2025.1_gpu, 2025.1
 dualsphysics          5.2.1_gpu, 5.2.1
 fds                   6.9.1, 6.8
 fvcom                 5.1.0
 gromacs               2025.0_gpu, 2025.0, 2022.2_gpu, 2022.2
 gx                    11-2024_gpu
 mohid                 24.10
 nwchem                7.2.3, 7.2.2
 openfast              4.0.2, 3.5.2
 openfoam-esi          2406, 2412, 2206
 openfoam-foundation   12, 8
 opensees              2.5.0, 3.7.1
 openseespy            3.7.1
 quantum-espresso      7.4.1, 7.3.1
 reef3d                25.02, 24.02, 24.12
 schism                5.11.0
 snl-swan              2.2
 splishsplash          2.13.0
 swan                  41.45, 41.31, 41.51
 swash                 10.05, 11.01, 10.01A, 9.01A
 xbeach                1.24, 1.23
```

## Flags
### `-h, --help`

Show help message and exit.

## Need Help?
Run the following command for more details:

```sh
inductiva simulators --help
```
