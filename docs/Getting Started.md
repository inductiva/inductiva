# Getting Started
Welcome to Inductiva, where large-scale simulations meet simplicity, through a 
seamless, user-friendly interface that brings cutting-edge computational resources 
to your fingertips—effortlessly.

This section is designed to provide you with all the necessary information you 
need to get started with Inductiva API. Here, we've condensed everything so you 
can install Inductiva, grab usage tips, and find some of its essential features, 
with ease.

It's all you need to get familiar and make the most out of Inductiva from the 
get-go.

In this section you will find:

* [Overview]()
    * [What is Inductiva?]()
    * [Who is it for?]()
    * [What can you do?]()
    * [Why choose Inductiva?]()
* [Prerequisites]()
* [Installing Inductiva]()
* [What to read next]()
* [Troubleshooting]()
* [Frequently Asked Questions (FAQs)]()
* [Useful Terminology]()


## Overview
### What is Inductiva?
Inductiva is a state-of-the-art Python package that simplifies running large-scale 
physical system simulations. Just two lines of code unleash its potentional, freeing 
you from local hardware limits and the technical hassles of scaling simulations. 
It's all about making big computational tasks a breeze.

### Who is it for?
Inductiva is perfect for curious minds in research, science, and engineering who 
are engaged in intensive computational work. Whether you're unraveling the mysteries 
of molecular dynamics, exploring the complexities of fluid dynamics, or solving 
complex engineering problems, Inductiva provides the computational power and simplicity 
you need. It's perfect for those who want to push the boundaries of simulation-based 
research without getting bogged down in the technicalities of simulation setup and 
execution. It turns brain-busting calculations into child's play, so you can spend 
more time having those 'eureka!' moments and less time waiting for your simulations 
to run.

### What can you do with Inductiva API?
With Inductiva, you have the flexibility to scale your simulations according to 
your needs:

* For a quick demo, there are [ready-to-use simulators]() that allow you to start 
simulations immediately with your existing configuration files. No fuss over 
installations or tangled software wires. 

* For those who dream big and simulate bigger, there's the option to [build your own]() 
digital dream team of machines, known as a ['MachineGroup'](). This feature lets 
you configure a custom-tailored squad of virtual machines ready to tackle those 
heavyweight simulations. 

So, whether you're in for a quick demo dash or for comprehensive, detailed analyses, 
Inductiva is equipped to facilitate your work seamlessly.

### Why choose Inductiva API?
Because choosing Inductiva means opting for simplicity and efficiency in your 
simulation work. With just two lines of code, you can activate its powerful features. 
No need to tangle with technical headaches or a web of dependencies. Inductiva 
streamlines the process, so can you cut through the clutter and focus on what 
really matters – *the science*.

## Prerequisites

### Python Installation
Inductiva requires Python version 3.8 or higher to function properly. If Python 
is not already installed on your system, follow these steps:

1. Visit the [Python download page](https://www.python.org/downloads/).
2. Select and download the appropriate Python version for your operating system.
3. Follow the installation instructions to set up Python on your system.

> **Tip:** Python updates are released annually in October. It's worth noting that 
it may take a few months for the scientific Python ecosystem to fully adapt to new 
releases. Therefore, for the best compatibility with Inductiva, consider using a 
Python version that has been stable for a few months.

### Visualization Tools
While using the [Inductiva API](), keep in mind that it does not provide built-in 
visualization tools for the simulators. To visualize your simulation results:

1. Download the results from Inductiva.
2. Utilize any visualization tool of your choice to analyze and interpret the data. 
For large-scale, complex 3D visualizations, our team recommends open-source tools 
like [ParaView](https://www.paraview.org/download/) or [VTK](https://vtk.org/download/).

## Installing Inductiva API

To install the latest Inductiva package release on PyPI with pip, just open your 
terminal and run:

```
pip install --upgrade inductiva
```

This will provide the core functionalities of the API, which is designed to facilitate 
various tasks in a computational environment. It will allow you to easily dispatch 
computational tasks, manage system resources, and run detailed simulations.

Encountering issues? Don’t worry. Head over to [Troubleshooting]() to work it out.

Otherwise, you're all set! You can now decide [what to read next]().

## What to read next

You can start experimenting with Inductiva by choosing from a variety of ready-to-use, 
[open-source simulators]() that are immediately available for use on our default 
machines.

You can set up custom, dedicated machine groups for more complex or [larger-scale simulation]() 
tasks, tapping into the more extensive [computational resources]() of Inductiva.

You can get inspired by how other researchers have used Inductiva API to configure, 
execute, and scale their simulations with minimal coding effort by visiting our 
[example gallery]().

You can troubleshoot installation problems that you might encounter with Inductiva 
API by checking [troubleshooting](#troubleshooting), [FAQs](), or [getting in touch]() 
with our support team directly.

## Troubleshooting

If installing the package fails, you can retry it on a new Python virtual environment. 
A [virtual environment](https://docs.python.org/3/library/venv.html) allows you to 
have a fresh Python environment with isolated dependencies. In your shell, run:

```
python -m venv <venv>
```

In that command, you should replace `<venv>` with the path (*e.g.*, `.venv`) in 
which you would like to create the environment. Then, to activate the environment 
(again, correctly replacing `<venv>`), run:

For `bash`/`zsh`:

```
source <venv>/bin/activate
```

For `cmd.exe` (Windows):

```
<venv>\Scripts\activate.bat
```

For `PowerShell` (Windows):
```
<venv>\Scripts\Activate.ps1
```

After activating the virtual environment, you can install the package as described 
below:

```
pip install --upgrade inductiva
```

## Frequently Asked Questions
This section could be removed or updated later.
## Useful Terminology
This documentation uses the following terms:

https://packaging.python.org/en/latest/glossary/ 
