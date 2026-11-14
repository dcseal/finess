FINESS plotting routines based on the MATrix LABoratory (MATLAB)
================================================================

This library contains plotting routines for [FINESS](../../README.md) that
are written in [MATLAB](http://www.mathworks.com/).  Many of these work
may work out of the box with the open-source alternative to MATLAB called 
[Octave](http://www.gnu.org/software/octave/), but the developers do not keep
up to date with this option.

See also: ['$FINESS/viz/python'](../python/README.md) for a complete set of 
[Python](http://www.python.org/)
plotting routines based on [matplotlib](http://matplotlib.org/).

NOTE: these routines were adapted from those found in the
[DoGPack](http://www.dogpack-code.org) software package.

Basic overview
--------------

From an application folder, after running `finess.exe`, call one of

    plotfin1(), plotfin2()

for 1D and 2D applications.  The basic format is

    plotfin[n]( parameters_ini_filename )

where `parameters_ini_filename` (default: `'parameters.ini'`) is the
parameters file used for the run.  The output directory is taken from
`output_dir` in the `[finess]` section of that file (default: `'output'`).

User supplied script
--------------------

Each application may supply a single script that is called once for each
frame that is plotted:

    plotq1 - 1D script
    plotq2 - 2D script

These routines allow the user to supply items such as line colors, axes, grid
labels, color bars, etc., everything that is application specific to the given
problem.

In the case of a user not supplying one of these function, a single generic
routine is supplied here in the library.
