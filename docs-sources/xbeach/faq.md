**Find answers to commonly asked questions about XBeach.**

<br>

# FAQ

## 1. Why does my XBeach task fail with the error: *Number of mpi domains in M-direction is greater than nx/4*?
XBeach requires that the number of MPI partitions in the M-direction (usually aligned with the x or nx 
dimension of your grid) does not exceed `nx / 4`. If you run your simulation on a machine with a large 
number of vCPUs (for example, a machine with 96 vCPUs), XBeach may attempt to create more partitions 
than the domain can accommodate, resulting in the error message:

```
Number of mpi domains in M-direction is greater than nx/4
XBeach cannot split into separate model domains.
```

For example, with a computational grid where:
- `nx = 360`
- `ny = 180`

The maximum number of MPI partitions in the M-direction is:
- `360 / 4 = 90`

If you use a machine with more than 90 vCPUs (e.g., `c4-highmem-96` with 96 vCPUs), XBeach will try to split the domain into more partitions than allowed, leading to the error.

### Why doesn't this happen on my laptop?
On a typical laptop, the number of available computational cores is much smaller than on high-performance machines, making it less likely to hit the partition limit. However, on larger machines or instances with more vCPUs, XBeach is more likely to run into this issue.

### Recommended Solution
To avoid this error, use a VM instance with fewer vCPUs. Ensure that the number of vCPUs you allocate does not exceed `nx / 4`. For instance, with `nx = 360`, limit the number of vCPUs to 90 or fewer.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
