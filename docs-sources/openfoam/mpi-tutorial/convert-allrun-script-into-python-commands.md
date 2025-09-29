# Convert `Allrun` into Python Commands for Multi-Node Execution
With Inductiva, running your simulations on a multi-node MPI (Message Passing Interface) cluster is as straightforward as 
using any single-node option. However, before proceding, there is one important step: you need to convert your `Allrun` 
script into Python commands.

To fully leverage MPI clusters, you **cannot** use the `shell_script` argument as you might with OpenFOAM simulations. 
Instead, you must use the `commands` argument, which accepts a list of commands to be executed during the simulation. 
Some commands will run sequentially, while others may run in parallel.

Converting your `Allrun` script into a list of commands is simple: extract each command from the script, excluding any logic or flow control, and place them as individual strings inside a Python list.

To illustrate, consider this shell snippet:

```bash
echo -e -n '\n    Running snappyHexMesh ...\n'
(
    cd system
    cp controlDict.SHM controlDict
    cp fvSchemes.SHM fvSchemes
    cp fvSolution.SHM fvSolution
    cd ..
)
```

This would become the following in Python:

```python
commands = [
    # You can also add the echo if you wanted
    "cp system/controlDict.SHM system/controlDict",
    "cp system/fvSchemes.SHM system/fvSchemes",
    "cp system/fvSolution.SHM system/fvSolution",
]
```

> ⚠️ **Important Notes**
>
> * Each command must be independent. You **cannot** use `cd` to change directories between commands.
> * Remove any logic such as conditionals or loops.
> * Special characters like `>`, `|`, and `&` are **not allowed**.