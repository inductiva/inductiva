# SWAN

SWAN is a simulator for obtaining realistic estimates of wave parameters in coastal
areas, lakes and eastuaries from given wind, sea floor and current conditions.

The simulator is configured using a single file with the `.swn` extension, and
additional files containing information about the domain and the sea floor, and
the conditions shall be saved in an input directory that is passed to the simulator.

## Example

```{literalinclude} ../../examples/swan/swan.py
:language: python
```

Check the [official documentation](https://swanmodel.sourceforge.io/) of SWAN to know 
more about the configuration details specific of the simulator.

## Inductiva Benchmarks

The following benchmark is currently available for SWAN:

* [Ring](https://benchmarks.inductiva.ai/SWAN/ring/): The "Ring" example from 
SWAN's site.

## What to read next

If you are interested in SWAN, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [Reef3D](Reef3D.md)
* [SCHISM](SCHISM.md)
* [SWASH](SWASH.md)
* [XBeach](XBeach.md)
