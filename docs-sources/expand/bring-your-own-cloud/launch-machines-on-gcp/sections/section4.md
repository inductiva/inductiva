# Troubleshooting and Limitations

## Current Limitations (Alpha Version)

Please be aware of the following limitations in the current alpha release:

### Machine Scaling
- **Single machine per group**: BYOC machine groups only support one VM per group (unlike regular machine groups where you can specify `num_machines`)
- **Workaround**: Create multiple machine groups in a loop to scale up (see [Python Client Usage](section2.md#machine-scaling-workaround))

### Cost Calculation
- **Cost reporting**: Currently always reports 0 in cost calculations

### Storage
- **Storage support**: Only Inductiva storage is currently supported
- **GCP storage integration**: Coming soon in future releases

## Troubleshooting

### Common Issues and Solutions

**Authentication Errors**
- Run `gcloud init` to set up authentication and project configuration
- Verify with `gcloud auth list` that you're authenticated
- Check your default project with `gcloud config get-value project`

**Permission Denied Errors**
- Ensure your Google account has the required permissions:
  - `compute.instances.create`
  - `compute.instances.delete` 
  - `compute.instances.setMetadata`
- Check your project's IAM settings in the [GCP Console](https://console.cloud.google.com/iam-admin/iam)

**Quota Exceeded Errors**
- Check your GCP project quotas in the [Quotas page](https://console.cloud.google.com/iam-admin/quotas)
- Request quota increases for the specific machine types you need
- Consider using different machine types or regions with available capacity

**Billing Issues**
- Ensure your GCP project has a valid billing account attached

**API Not Enabled**
- Enable the Compute Engine API: `gcloud services enable compute.googleapis.com`

**Machine Not Starting**
- Check the GCP Console for VM instance status
- Verify that the task-runner container is running
- Check VM logs in the GCP Console for any startup errors

**Network Connectivity Issues**
- Ensure your GCP project allows outbound HTTPS traffic (port 443)
- Check if your organization has firewall rules blocking external connections
- Verify that the VM has internet access


```{banner_small}
:origin: launch_machines_on_gcp_sec4
```
