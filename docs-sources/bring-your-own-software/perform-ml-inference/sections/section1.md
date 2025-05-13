# DiffDock Setup 🔧
[DiffDock](https://github.com/gcorso/DiffDock) is an innovative molecular docking tool that uses diffusion models to predict protein-ligand interactions. 
Unlike traditional docking methods, DiffDock uses a generative AI approach to predict the binding of small molecules to proteins with high accuracy and speed. 
This model has shown great promise for drug discovery and virtual high throughput screening.

> For a deeper understanding, explore the original paper 👉 ["DiffDock: Diffusion steps, twists and turns for molecular docking"](https://arxiv.org/abs/2210.01776).

A key challenge with DiffDock is that the default `rbgcsail/diffdock:latest` Docker image pulls model checkpoints from Hugging Face or other remote sources at runtime. 
This causes significant delays when running inferences, as model downloads occur each time a container is spun up. 

To resolve this, we can pre-download the models and bundle them into a custom Docker image. This approach offers several advantages:
- **Eliminate runtime model downloads**
- **Improve inference startup time**
- **Make the image self-contained and portable**

Let’s walk through the process step by step.

## Step 1: Identify Downloaded Models 
Before building a custom image, we need to determine:
- Which models are being downloaded
- Where they're stored in the container

### How to Inspect Model Downloads 
1. Run a test container using the official image:
```bash
docker run --rm -it rbgcsail/diffdock:latest bash
```

2. Run the inference script (or partial commands) to trigger model downloads.

3. Check the common model cache paths:
```bash
ls -lh ~/.cache/huggingface
ls -lh ~/.torch/hub/checkpoints
```

In the case of DiffDock, the models typically reside in the following locations:
`./hub/checkpoints/` (for ESM models from Hugging Face)
`./workdir/v1.1/score_model/`
`./workdir/v1.1/confidence_model/`

Take note of these files, as they will be manually downloaded and included in the custom image.

### Model Files to Include
From a test run, we identified the following key models to cache:

**Confidence Model** (`confidence_model`):
- `best_model_epoch75.pt`
- `model_parameters.yml`

**Score Model** (`score_model`):
- `best_ema_inference_epoch_model.pt`
- `model_parameters.yml`

**ESM Models** (in hub/checkpoints):
- `esm2_t33_650M_UR50D.pt`
- `esm2_t33_650M_UR50D-contact-regression.pt`
- `esm2_t36_3B_UR50D.pt`
- `esm2_t36_3B_UR50D-contact-regression.pt`
- `esmfold_3B_v1.pt`

You will organize these models locally as follows:

```bash
.
├── confidence_model/
├── score_model/
└── esm_models/
```

These will later be copied into the container.

## Step 2: Write a Custom Dockerfile
In this step, we'll create a Dockerfile that:

- Starts from the official `rbgcsail/diffdock:latest` image
- Copies your pre-downloaded models into their correct locations

As shown below, here's an example of the Dockerfile:

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
ENV TORCH_HOME=/home/appuser/

CMD ["bash"]
```

## Step 3: Build the Docker Image
Once your models and Docker file are in the same directory, build the image:

```
docker build -t diffdock-with-models .
```

This can take a few minutes, depending on your system and the size of your models. 
Feel free to grab a coffee while you wait! ☕
