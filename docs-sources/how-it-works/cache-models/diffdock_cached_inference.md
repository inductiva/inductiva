# 🚀 Accelerating DiffDock Inference with a Custom Docker Image on Inductiva

This guide outlines how to create a custom Docker image with preloaded machine learning models, ready to run on the Inductiva API. The goal is to reduce model load time during inference by caching all necessary assets inside the container image.

We use the [DiffDock](https://github.com/gcorso/DiffDock) model as a working example.

--- 

## 🎯 Goal
The default rbgcsail/diffdock:latest Docker image pulls model checkpoints from Hugging Face or other remote sources at runtime. This causes significant delays when running inferences, as model downloads happen each time a container is spun up.

By pre-downloading the models and bundling them into a custom image, we:

- Eliminate runtime model downloads
- Improve inference startup time
- Make the image self-contained and portable

--- 

## 🔍 Step 1: Identify Downloaded Models
Before building a custom image, we need to determine:

- Which models are being downloaded
- Where they're stored in the container


### 🔎 How to inspect model downloads:

1. Run a test container with the official image:

```bash
docker run --rm -it rbgcsail/diffdock:latest bash
```

2. Run the inference script (or partial commands) to trigger model downloads.

3. Look inside the common model cache paths:

```bash
ls -lh ~/.cache/huggingface
ls -lh ~/.torch/hub/checkpoints
```

In DiffDock’s case, the models typically go to:
- ./hub/checkpoints/ (used for ESM models from Hugging Face)
- ./workdir/v1.1/score_model/
- ./workdir/v1.1/confidence_model/

Take note of these files — they will be manually downloaded and included in the custom image.


--- 

## 📁 Model Files to Include

From a test run, we identified the following key models to cache:

### Confidence Model (confidence_model)

- best_model_epoch75.pt
- model_parameters.yml

### Score Model (score_model)

- best_ema_inference_epoch_model.pt
- model_parameters.yml

### ESM Models (in hub/checkpoints)

- esm2_t33_650M_UR50D.pt
- esm2_t33_650M_UR50D-contact-regression.pt
- esm2_t36_3B_UR50D.pt
- esm2_t36_3B_UR50D-contact-regression.pt
- esmfold_3B_v1.pt

You’ll organize these locally into:

```bash
.
├── confidence_model/
├── score_model/
└── esm_models/
```

These will be copied into the container later.

---

## 🐳 Step 2: Write a Custom Dockerfile

We’ll create a Dockerfile that:
- Starts from the official rbgcsail/diffdock:latest image
- Copies your pre-downloaded models into their correct locations

### 🧾 Dockerfile

```Dockerfile
# Use the base DiffDock image
FROM rbgcsail/diffdock:latest

WORKDIR /home/appuser/

# Create model directories
RUN mkdir -p ./workdir/v1.1/confidence_model \
             ./workdir/v1.1/score_model \
             ./hub/checkpoints

# Copy model checkpoints and assign ownership directly
COPY --chown=appuser:appuser confidence_model/ ./workdir/v1.1/confidence_model/
COPY --chown=appuser:appuser score_model/ ./workdir/v1.1/score_model/
COPY --chown=appuser:appuser esm_models/ ./hub/checkpoints/

# Set environment variable to avoid external downloads
# The TORCH_HOME will look for the hub/checkpoints directory
# and the ESM models will be found there
ENV TORCH_HOME=./

CMD ["bash"]
```

---

## 🔨 Step 3: Build the Docker Image
Once your models and Dockerfile are ready in the same directory, build the image:

```
docker build -t diffdock-with-models .
```

This may take a few minutes, depending on your system and the size of the models.

---


## 📦 Step 4: Export the Docker Image to Inductiva Storage

Once your Docker image is built, the next step is to upload it to your Inductiva storage so it can be used in remote simulations.



### 🧬 Upload the Image

The upload process is streamlined using the Inductiva CLI. One command handles everything — including conversion to the required .sif format:


```bash
inductiva containers upload diffdock-with-models
```

This command will:
- Locate the Docker image tagged diffdock-with-models
- Convert it to the Singularity Image Format (SIF)
- Upload it to your Inductiva containers registry

✅ Note: This requires the Inductiva CLI to be installed and properly configured. If you haven't set it up yet, follow the Inductiva install guide.

You can read more about the personal docker usage in the [Inductiva documentation](https://inductiva.ai/guides/documentation/intro/private-dockers).

### 📂 Verify the Upload

After the image has been uploaded, you can verify it's available in your account by listing your containers:

```bash
inductiva containers list
```

You should see diffdock-with-models listed among your available container images.


---

## 🚀 Step 5: Run the Custom Image on the Inductiva API

With your custom container uploaded to Inductiva, you're now ready to use it for simulations — no more model downloads at runtime 🎉.

### 🧪 Using the Image in a Simulation

You can specify your uploaded container in an `inductiva.simulators.CustomImage` object. Here’s an example of how to define and run a task using your `diffdock-with-models` image:

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
        'micromamba run -n diffdock python ./DiffDock/inference.py '
        '--config ./custom_model_path_inference_args.yaml '
        '--protein_sequence GIQSYCTPPYSVLQDPPQPVV '
        '--ligand "COc(cc1)ccc1C#N" '
        '--samples_per_complex 5'
    ]
)
```

> 💡 Make sure any required config files or inputs (like custom_model_path_inference_args.yaml) 
> are included in your input directory.


## ✅ Wrapping Up

With this workflow, you've created a fast, portable, and scalable version of DiffDock for use with Inductiva. By preloading models in a custom Docker image and uploading it to your personal container registry, you’ve eliminated costly download steps during inference — leading to faster, more reliable cloud simulations.

You’re now equipped to use custom deep learning workloads with maximum efficiency on Inductiva 🚀.