# Export the Docker Image to Inductiva Storage 📦
After building your Docker image, the next step is to upload it to your Inductiva storage, making it ready for use in remote simulations.

## Upload the Image
The upload process is streamlined by the Inductiva CLI. This *one command* does it all — including converting to the required `.sif' format:

```bash
inductiva containers upload diffdock-with-models
```

This command will:
- Find the Docker image tagged `diffdock-with-models`
- Convert it to the Singularity Image Format (SIF)
- Upload it to your Inductiva containers registry

> **Note**: Ensure the Inductiva CLI is installed and properly configured. If you haven't done so yet, refer to the Inductiva installation guide.

For more information on using personal Docker images, visit the [Inductiva documentation](https://inductiva.ai/guides/documentation/intro/private-dockers).

## Verify the Upload
Once the image is uploaded, you can confirm it's available in your account by listing your containers:

```bash
inductiva containers list
```

You should see `diffdock-with-models` listed among your available container images.

