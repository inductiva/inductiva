---
orphan: true
---

# Downloading The Results

We have successfully executed 50 simulations in parallel, each corresponding to
a different water depth variation. This demonstrates the power of leveraging
cloud resources to scale up computational experiments efficiently. Now, it's
time to retrieve all the results and analyze the data to extract meaningful
insights from our simulations.

| Simulation with water depth of 100 | Simulation with water depth of 200 |
|:---------------:|:-----------------:|
| <img src="../_static/openfast_animation_30_fps_100.gif" alt="OpenFAST simulation visualization"> |<img src="../_static/openfast_animation_30_fps_200.gif" alt="OpenFAST simulation visualization">|

## Project summary and download outputs

Using the Inductiva package, obtaining a project summary and downloading all
output files is straightforward. The following Python snippet demonstrates
how to achieve this:

```python
import inductiva

openfast_project = inductiva.projects.Project(
    name="Openfast_WavesWN")

print(openfast_project)

openfast_project.download_outputs()
```

Executing `print(openfast_project)` provides a summary with key project details:

```
Project 'Openfast_WavesWN' with 50 tasks (id=3dde0355-0189-4634-83eb-8e177e873bf3).

Tasks status:
success: 50

Total number of output files: 4150
Total size of output: 1.63 GB

Project duration: 7 minutes and 4 seconds
Project total simulated time: 28 minutes and 26 seconds

Estimated project cost: 0.0059 US$
```

### Key Takeaways

- **Efficiency Gains**: Running all 50 simulations in parallel took just
**7 minutes and 4 seconds**, compared to **28 minutes and 26 seconds** if
executed sequentially (~34.12 seconds per simulation).
- **Comprehensive Output**: The simulations generated **4,150 files** with a
total size of **1.63 GB**.
- **Cost-Effectiveness**: The entire project execution cost only **$0.0059**,
showcasing the affordability of leveraging cloud computing for large-scale simulations.

By efficiently managing computational resources, we achieve faster results while
keeping costs remarkably low.

Lastly, running the `openfast_project.download_outputs()` will create a folder
called `inductiva_output` with one folder for each simulation, as seen here:

```
total 0
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 1d61519z29drg2jp4nhknauoo
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 3qoor8iooxtmmm5dozs0o40px
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 4jho5mrt1dvr1c4wcsq9hje2z
...
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 tif6v2qugbkcjcum2hna5ttow
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 tyqb7rnsyzi94wqgcn58cc2y6
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 u3my9tc6deh72xan9j7kn0syj
```

In conclusion, leveraging cloud computing for large-scale simulations not only
enhances efficiency but also significantly reduces computational time and costs.
By running 50 simulations in parallel, we managed to work around one of OpenFAST
limitations, witch is the fact that it does not scale with multiples CPU cores.

Inductiva makes it possible and convenient to run hundreds thousands of
simulations. For example, you could now change the code for:

1. Covering a design space with many more parameters;
2. Adding a black box optimizer, such as [Vizier](https://github.com/google/vizier),
on top of the loop and adding a evaluation function to analyze the result of
each simulation that finishes, and perform an intelligent exploration of the
design space.
3. Building a really large dataset of example simulations to later train a
surrogate model;

Inductiva can streamline research, making high-performance computing more
accessible and cost-effective.
