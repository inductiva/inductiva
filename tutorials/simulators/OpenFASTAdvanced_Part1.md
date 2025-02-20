# OpenFAST Tutorial (Advanced)

In this tutorial, we are going to show you how you can
use the Inductiva API to accelerate your OpenFast projects.

## The peculiar case of OpenFast

OpenFast is a pretty peculiar piece of software because
it runs only on a single thread. That means that OpenFast
does not take full advantage of the extra compute capacity
offered by modern CPUs, which support dozens, or sometimes 
hundreds, of computational threads in parallel. In fact,
the most decisive factor impacting the performance of 
OpenFast simulations is CPU clock frequency. 

Now, what is interesting about that is that, tipically, your 
500\$ home desktop has a much higher clock frequency than the 
20000\$ monster machines available in cloud centers. Your local
desktop has CPU with clock frequencies between 4 and 5GHz (but
can be overclocked to more than 5.5GHz), while cloud machines
have CPUs that will be running at around 3GHz (or less!).

But why is that? Have you noticed how fast the ventilation 
fans of your desktop PC spin when you are running some long 
simulation? Now imagine how many fans and sophisticated
cooling systems you would need if you packed 100 CPUs like
yours in a volume the size of your kitchen fridge. Hardware
density in datacenter cabinets is so high that if all the
CPUs were running at speeds of 5GHz or more, the datacenter
would be impossible to cool down.

So, if that is the case, why would you even bother using
Inductiva to run your OpenFast simulation? What advantage
could you possibly have if you decide to send your OpenFast
simulation to a cloud machine using Inductiva? 

And the answer is: **hundreds of advantages!**

Let me explain! If you are running only 1 OpenFast simulation, 
then you should run it on your desktop machine: it will be faster
there due to much higher clock speeds.  But if you need to run
hundreds or thousands of OpenFast simulations, then you can use
Inductiva to spin up hundreds of very cheap cloud machines to 
run all those simulations in parallel instead of having to run
them sequentially in your machine! This will be hundreds of times
faster. And it is super easy (and cost-effective) with Inductiva

<h3>🌊 OpenFAST Simulation Visualization</h3>
<p align="center">
  <img src="../_static/openfast_animation_30_fps.gif" alt="OpenFAST simulation visualization" width="700">
</p>

In the pages to come you can see the following topics:
- asdasd
- asdasdasd-
- sadasda

[Preparation](OpenFASTAdvanced_Part2.md)