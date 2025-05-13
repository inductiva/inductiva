# Run the Custom Image on the Inductiva API ⚡️
With your custom container uploaded to Inductiva, you're now ready to run simulations — no more waiting for model downloads at runtime! 🎉

## Using the Image in a Simulation
You can specify your uploaded container in an `inductiva.simulators.CustomImage` object. Here's an example of how to define and run a task using your `diffdock-with-models` image:

```python
import inductiva

# Define the machine group
mg = inductiva.resources.MachineGroup(
    machine_type="g2-standard-8",
    num_machines=1,
    spot=True,
)

# Define the simulator using your custom image
diffdock_simulator = inductiva.simulators.CustomImage(
    container_image="inductiva://my-containers/diffdock-with-models"
)

# Define input parameters and run the task
task = diffdock_simulator.run(
    input_dir="path/to/input",
    on=mg,
    commands=[
        'micromamba run -n diffdock python /home/appuser/DiffDock/inference.py '
        '--config ./custom_model_path_inference_args.yaml '
        '--protein_sequence GIQSYCTPPYSVLQDPPQPVV '
        '--ligand "COc(cc1)ccc1C#N" '
        '--samples_per_complex 5'
    ]
)
```

> 💡 **Tip**: Ensure that any necessary config files or inputs (like `custom_model_path_inference_args.yaml`)
> are included in your input directory.

### Experiments on samples per complex

In order to explore the effect of varying the `samples_per_complex` parameter on runtime and cost, we conducted a series of benchmark tests using two machine types: `a2-highhpu-1g` (NVIDIA A100) and `g2-standard-8` (NVIDIA L4). This parameter controls the number of binding poses generated per protein-ligand pair, which directly impacts computational load.

The following table summarizes the duration and cost for processing different numbers of samples on each machine type:


| Machine Type    | Samples | Duration    | Cost   | Samples/sec | Cost/sample |
|-----------------|---------|-------------|--------|-------------|-------------|
| a2-highhpu-1g    | 100     | 17 min 15 s | $0.43  | 0.0966      | 0.004300    |
| a2-highhpu-1g    | 1000    | 22 min 25 s | $0.55  | 0.7435      | 0.000550    |
| a2-highhpu-1g    | 5000    | 45 min 35 s | $1.13  | 1.8282      | 0.000226    |
| g2-standard-8    | 100     | 17 min 39 s | $0.10  | 0.0944      | 0.001000    |
| g2-standard-8    | 1000    | 21 min 56 s | $0.12  | 0.7599      | 0.000120    |
| g2-standard-8    | 5000    | 46 min 18 s | $0.26  | 1.8008      | 0.000052    |


These results indicate that while `a2-highhpu-1g` achieves marginally faster runtimes, the `g2-standard-8` offers substantially better cost-efficiency, especially at larger sample sizes. Therefore, for many use cases, `g2-standard-8` may be the preferred option when optimizing for price-performance trade-offs.

## Wrapping Up
With this workflow, you've created a fast, portable, and scalable version of DiffDock for use with Inductiva. By preloading models in a custom Docker image and uploading it to your personal container registry, you've eliminated costly download steps during inference, resulting in faster and more reliable cloud simulations.

You're now equipped to run custom deep learning workloads with maximum efficiency on Inductiva! ✅

