# Uninstall the `inductiva` package

Whether you just want a completely new installation of `inductiva` or
simply want to completely remove the package from your system, we've
got you covered. Either way, we would like to give you instructions on
how to completely remove any trace of `inductiva` from your system.

⚠️**NOTE** No simulation files generated through the package will be
removed from your cloud storage. If you are not going to use Inductiva 
anymore, please delete your data. 

First, you should uninstall the package using:

```
pip uninstall inductiva
```

Next, recall that the inductiva package, writes locally to your
disk. Therefore, you should have a folder under `$HOME/.inductiva` in
your system. Remove it and you are good to go!

```
rm -rf $HOME/.inductiva
```