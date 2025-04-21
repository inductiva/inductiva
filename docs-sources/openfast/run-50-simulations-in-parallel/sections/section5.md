# Results and Key Takeaways
We successfully ran 50 simulations in parallel, each corresponding to a different water depth variation. This demonstrates the power of using cloud resources 
to efficiently scale up computational experiments. Now it's time to retrieve all the results and analyze the data to extract meaningful insights from our simulations.

| Simulation with water depth of 100 | Simulation with water depth of 200 |
|:---------------:|:-----------------:|
| <img src="./_static/openfast_animation_30_fps_100.gif" alt="OpenFAST simulation visualization"> |<img src="./_static/openfast_animation_30_fps_200.gif" alt="OpenFAST simulation visualization">|

## Project Summary and Output Download
Using the Inductiva package, it is easy to get a project summary and download all the output files. The following code snippet demonstrates how to do this:

```python
import inductiva


openfast_project = inductiva.projects.Project(
   name="Openfast_WavesWN")


print(openfast_project)


openfast_project.download_outputs()
```


Executing `print(openfast_project)` gives a summary of the main project details:


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

## Key Takeaways
- **Efficiency gains**: Running all 50 simulations in parallel took only 7 minutes and 4 seconds, compared to 28 minutes and 26 seconds when
run sequentially (~34.12 seconds per simulation).
- **Comprehensive output**: The simulations generated 4,150 files with a total size of 1.63 GB.
- **Cost-effectiveness**: The entire project cost only $0.0059 to run, demonstrating the affordability of using cloud computing for
large-scale simulations. By efficiently managing computing resources, we achieve faster results while keeping costs remarkably low.

Finally, running `openfast_project.download_outputs()` creates a folder called `inductiva_output/Openfast_WavesWN` with one folder for each simulation, as shown here:

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

In summary, using cloud computing for large-scale simulations not only increases efficiency, but also significantly reduces computational time and cost. 
By running 50 simulations in parallel, we were able to overcome one of the limitations of OpenFAST, which is that it does not scale with multiple CPU cores.

Inductiva makes it possible and convenient to run hundreds or thousands of simulations. For example, you could now change the code for:

1. Cover a design space with many more parameters;
2. Add a black-box optimizer, such as [Vizier](https://github.com/google/vizier),
on top of the loop and add an evaluation function to analyse the result of each completed simulation and perform an intelligent exploration of the design space;
3. Build a really large dataset of example simulations to later train a surrogate model.

Inductiva can simplify research by making high-performance computing more accessible and cost-effective.