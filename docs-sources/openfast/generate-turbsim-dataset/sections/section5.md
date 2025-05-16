# Results and Key Takeaways
We successfully generated a TurbSIM dataset in parallel by sampling the seeds and the wind speed.
This demonstrates the power of using cloud resources to efficiently scale up computational experiments. 
Now it's time to retrieve all the results and analyze the data to extract meaningful insights from our simulations.

<p align="center"><img src="../../_static/turbsim_animation_30_fps.gif" alt="TurbSIM simulation visualization" width="700"></p>


## Project Summary and Output Download
Using the Inductiva package, it is easy to get a project summary and download all the output files. The following code snippet demonstrates how to do this:

```python
import inductiva


turbsim_project = inductiva.projects.Project(
   name="turbsim_dataset")


print(turbsim_project)


turbsim_project.download_outputs()
```


Executing `print(turbsim_project)` gives a summary of the main project details.

Finally, running `openfast_project.download_outputs()` creates a folder called `inductiva_output/turbsim_project` with one folder for each simulation.


## Key Takeaways
In summary, using cloud computing for large-scale simulations not only increases efficiency, but also significantly reduces computational time and cost. 
By running 50 simulations in parallel, we were able to overcome one of the limitations of TurbSIM, which is that it does not scale with multiple CPU cores.

Inductiva makes it possible and convenient to run hundreds or thousands of simulations. For example, you could now change the code for:

1. Cover a design space with many more parameters;
2. Add a black-box optimizer, such as [Vizier](https://github.com/google/vizier),
on top of the loop and add an evaluation function to analyse the result of each completed simulation and perform an intelligent exploration of the design space;
3. Build a really large dataset of example simulations to later train a surrogate model.

Inductiva can simplify research by making high-performance computing more accessible and cost-effective.
