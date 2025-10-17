**Find answers to commonly asked questions about WAVEWATCH III.**

<br>

# FAQ

## 1. What values can I use for the `switch` argument?
We provide all the default switches that WAVEWATCH III supports. So, you can use any of the following values for the `switch` argument:

|            |             |               |             |            |
|------------|-------------|---------------|-------------|------------|
| NCEP_gwm   | NRL2        | OASOCM        | UKMO_uk     | UoM_nl3s   |
| Ifremer1   | NCEP_st2    | NRL3          | SMCMlt      | USACE_1    |
| ite_pdlib  | Ifremer2    | NCEP_st4      | NRL4        | UKMO       |
| USACE_2    | multi_esmf  | Ifremer2_pdlib| NCEP_st4sbs | OASACM     |
| UKMO_gbl   | UoM_nl1     | swin          | NCEP_glwu   | NRL1       |
| OASICM     | UKMO_reg    | UoM_nl3       | ugdev2      |            |

You can also provide your custom switch in your input files by pointing to it using the `custom_switch` argument.

> **Note**: You can only use either the switch argument or the custom_switch argument — not both.

<br>

## 2. I keep getting the error `PDLIB requires METIS or SCOTCH library for domain decomposition`. How can I fix it?

This error usually occurs because the **switch file** used to compile the simulator is missing the required domain decomposition library.

For example, your switch file might contain something like:

```
PDLIB O2b
```

In this configuration, **PDLIB** is enabled without **METIS** or **SCOTCH**, which are required libraries for domain decomposition.

To fix the issue, simply update your switch file to include one of these libraries, for example:

```
PDLIB METIS O2b
```

This tells the compiler to build **PDLIB** with **METIS** support, resolving the error.


<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
