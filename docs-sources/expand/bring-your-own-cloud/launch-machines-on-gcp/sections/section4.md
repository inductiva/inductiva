# Troubleshooting and Limitations

## Current Limitations (Alpha Version)

Please be aware of the following limitations in the current alpha release:

### Machine Scaling
- **Single machine per group**: BYOC machine groups only support one VM per group (unlike regular machine groups where you can specify `num_vms`)
- **Workaround**: Create multiple machine groups in a loop to scale up (see [Python Client Usage](section2.md#machine-scaling-workaround))

### Cost Calculation
- **Cost reporting**: Currently always reports 0 in cost calculations

### Storage
- **Storage support**: Only Inductiva storage is currently supported
- **GCP storage integration**: Coming soon in future releases

## Troubleshooting

### Authentication Errors

**Problem**: Getting authentication errors when trying to launch machines.

**Solution**: Ensure your GCP credentials are properly configured:
- Run `gcloud init` to set up authentication and project configuration
- Verify with `gcloud auth list` that you're authenticated
- Check that your default project is set with `gcloud config get-value project`
- If using service accounts, ensure the service account key is properly configured

### Insufficient Permissions

**Problem**: Getting permission denied errors when creating VMs.

**Solution**: Ensure your Google account has the necessary IAM role:
- Compute Instance Admin (v1) role
- Check your project's IAM settings in the [GCP Console](https://console.cloud.google.com/iam-admin/iam)

### Insufficient Quotas

**Problem**: Getting quota exceeded errors when launching machines.

**Solution**: 
- Check your GCP project quotas in the [Quotas page](https://console.cloud.google.com/iam-admin/quotas)
- Request quota increases for the specific machine types you need
- Consider using different machine types or regions with available capacity

### Billing Issues

**Problem**: Getting billing-related errors.

**Solution**: Ensure your GCP project has a valid billing account attached.

### API Not Enabled

**Problem**: Getting errors about Compute Engine API not being available.

**Solution**: Enable the Compute Engine API:
```bash
gcloud services enable compute.googleapis.com
```

### Machine Not Starting

**Problem**: Machine group creation succeeds but the VM doesn't start properly.

**Solution**: 
- Check the GCP Console for VM instance status
- Verify that the task-runner container is running
- Check VM logs in the GCP Console for any startup errors
- Ensure your project has sufficient resources and quotas

### Network Connectivity Issues

**Problem**: Task-runner can't connect to Inductiva's backend.

**Solution**:
- Ensure your GCP project allows outbound HTTPS traffic (port 443)
- Check if your organization has firewall rules blocking external connections
- Verify that the VM has internet access

## Getting Help

If you continue to experience issues:

1. Check the [GCP Console](https://console.cloud.google.com/) for VM status and logs
2. Verify your gcloud CLI configuration with `gcloud config list`
3. Test basic GCP functionality with `gcloud compute instances list`
4. Contact Inductiva support with specific error messages and your GCP project details

```{banner_small}
:origin: launch_machines_on_gcp_sec4
```
