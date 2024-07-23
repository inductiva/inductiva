# NWChem

NWChem is a powerful and versatile computational chemistry software designed for
simulating molecular and materials systems at a high level of accuracy. Developed
by the Environmental Molecular Sciences Laboratory (EMSL) at Pacific Northwest
National Laboratory (PNNL), NWChem provides a wide range of quantum mechanical
methods, including density functional theory (DFT), Hartree-Fock, and various
post-Hartree-Fock approaches. It is capable of performing large-scale simulations
on modern high-performance computing architectures, making it suitable for
tackling complex chemical phenomena in both isolated molecules and extended
materials. NWChem's modular design allows for the seamless integration of new
functionalities, fostering ongoing development and collaboration within the
scientific community. This flexibility and robustness have made NWChem a go-to
tool for researchers seeking to explore chemical processes with high computational
demands.

## Example

```python
import inductiva

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "nwchem-input-example.zip", unzip=True)

nwchem = inductiva.simulators.NWChem()

task = nwchem.run(input_dir=input_dir,
                  sim_config_filename="h2o_sp_scf.nw",
                  n_vcpus=1)

task.wait()
task.download_outputs()
```
