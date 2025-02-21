---
orphan: true
---

# Downloading The Results

We have successfully executed 40 simulations in parallel, each corresponding to
a different water depth variation. This demonstrates the power of leveraging
cloud resources to scale up computational experiments efficiently. Now, it's
time to retrieve all the results and analyze the data to extract meaningful
insights from our simulations.

| |  |
|---------------|-----------------|
| <p align="center"><img src="../_static/openfast_animation_30_fps_180.gif" alt="OpenFAST simulation visualization" width="700"></p> |<p align="center"><img src="../_static/openfast_animation_30_fps_220.gif" alt="OpenFAST simulation visualization" width="700"></p>|

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
Project 'Openfast_WavesWN' with 40 tasks (id=3bb3f544-4d08-46d7-98de-6b1f3d55e5a9).

Tasks status:
success: 40

Total number of output files: 3320
Total size of output: 1.31 GB

Project duration: 5 minutes and 17 seconds
Project total simulated time: 24 minutes and 14 seconds

Estimated project cost: 0.0050 US$
```

### Key Takeaways

- **Efficiency Gains**: Running all 40 simulations in parallel took just
**5 minutes and 17 seconds**, compared to **24 minutes and 14 seconds** if
executed sequentially (~36.35 seconds per simulation).
- **Comprehensive Output**: The simulations generated **3,320 files** with a
total size of **1.31 GB**.
- **Cost-Effectiveness**: The entire project execution cost only **$0.0050**,
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
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 4rzzrir5dmhu4g49ft2t8dmw0
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 4vocm7rjiehyiyb05mmnpxw0u
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 5m9nnjw28h31puxybbzk1mq6h
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:10 60ygm64qsd56fmcjxgytq93pr
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 68ksaoricizknyhyg65r6cuqo
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 691acms5avqwvpcqf8kvecw33
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 7jkebn05rfsfcprto45k8nb4h
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 8gvt5i0x0qeaj2x7ftf7pwxyd
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:10 9ousr15npia4f9fphfjm1rmrn
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 b0mhrbzy9fwrfuj3f84q8ou75
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 bfn4ok5cwy93bmi1hg6vvk0rs
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 btctwrbp2z69db7w5rbdyuatm
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 bzhd49wxjkroc30f4vgs8u4o8
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 d6ght4mi7u10mzkstm33fpvty
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 eyw55mcvecqhjxnwncpy6wqkm
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 fr286c6cfrmcus12nh4n37m78
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 gnq3smc96ht404hz4ucvoovw6
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:10 h4qtfsq43si18a6oyvj74qrce
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 hmf133908834ms1hg9oh49q4e
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 iox5b1040tagfnwi9a3j26axr
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 j85bg918c1i0obj1izumtwwzw
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 jxn9a281jkyi4f62n4w80a1yp
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 k6xoup4jqczhjuwmiig4c7eo9
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 kndw77s3jzumypsnssszp4m4h
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 n0jb7r2rpuyik5xr7e96d024r
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 pqv4qmy2so820yoqm6dj51na0
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:10 pt3i70i4hmzn740og3ber5yvs
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 qfrv5vkims0zqbv7gevviiofg
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 qgehpoemdijghnu27hqjed813
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 r1kqwy76fvus9zaevv9xzzc02
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 smnphh002n0pxqnte7r6ihekg
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:10 snj7xr01k099qbh9ytf97qtbd
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:10 sy2a5r3cby0in6151ze06y29g
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 t98ihj29ieep26ndstns63u9r
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 tif6v2qugbkcjcum2hna5ttow
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 tyqb7rnsyzi94wqgcn58cc2y6
drwxr-xr-x  3 paulobarbosa  staff  96 Feb 20 15:11 u3my9tc6deh72xan9j7kn0syj
```

In conclusion, leveraging cloud computing for large-scale simulations not only
enhances efficiency but also significantly reduces computational time and costs.
By running 40 simulations in parallel, we managed to work around one of OpenFAST
limitations, witch is the fact that it does not scale with multiples CPU cores.
This was just a small example, and the same could be done to hundreds of simulations.

Inductiva can streamline research, making high-performance computing more accessible and cost-effective.
