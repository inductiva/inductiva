# Uninstall the `inductiva` package

Whether you just want a completely new installation of `inductiva` or
simply want to completely remove the package from your system, we've
got you covered. Either way, we would like to give you instructions on
how to completely remove any trace of `inductiva` from your system.

**NOTE** No simulation files generated though the package will be
removed.

First, you should uninstall the package using:

```shell
pip uninstall inductiva
```

Next, recall that the inductiva package, writes locally to your
disk. Therefore, you should have a folder under `$HOME/.inductiva` in
your system. Remove it and you are good to go!

```shell
rm -rf $HOME/.inductiva
```