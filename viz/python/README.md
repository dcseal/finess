FINESS plotting routines based on Python's MATPLOTLIB
======================================================

This library contains plotting routines for [FINESS](../../README.md) that
are written in [Python](https://www.python.org/) using
[matplotlib](https://matplotlib.org/).  They require Python 3 with numpy and
matplotlib installed.

See also: ['$FINESS/viz/matlab'](../matlab/README.md) for the
MATLAB plotting routines.

Basic overview
--------------

After running `finess.exe`, type

    $ python3 $FINESS/viz/python/plotfin[n].py

from the application directory, where [n] = 1 or 2 is the dimension of the
application.  The script asks which component of q to plot, then which frame;
press return to step to the next frame, or type `q` to quit.

For a full list of options, type

    $ python3 $FINESS/viz/python/plotfin[n].py -h

The options are

    -o, --outputdir   output directory to read from      (default: output)
    -q, --q-name      prefix of the solution files       (default: q)
    -x, --aux-name    prefix of the aux files            (default: a)

For example,

    $ python3 $FINESS/viz/python/plotfin1.py -o output_RK_9

plots the 1D solution stored in the directory `output_RK_9`.

User supplied script
--------------------

Each frame is drawn by a single function that an application may supply in
its own directory.  If the application directory contains this file, it is
used; otherwise the default in this directory is used.

| Dimension | Local file  | Function signature                                              | Default             |
|-----------|-------------|-----------------------------------------------------------------|---------------------|
| 1D        | `plotq1.py` | `plotq1(xc, qsoln, auxsoln, time, ini_params, m)`               | `plotq1_default.py` |
| 2D        | `plotq2.py` | `plotq2(xc, yc, xl, yl, qsoln, auxsoln, time, ini_params, m)`   | `plotq2_default.py` |

where

* `xc` (and `yc`) are the cell centers; in 2D, `xl`, `yl` are the cell
  corners with shape `(mx+1, my+1)`, suitable for `pcolor`.
* `qsoln` has shape `(mx, meqn)` in 1D and `(mx, my, meqn)` in 2D.
* `auxsoln` has the same layout with `maux` components, or is `0.0` when
  `maux = 0`.
* `ini_params` is a dictionary of the values in `output/parameters.ini`.
* `m` is the zero-based component of q that was selected.

These routines let the user supply items such as line colors, axes, grid
labels, color bars and exact solutions; everything that is specific to the
application.  See `$FINESS/apps/2d/advection/sine_example/plotq2.py` for a
2D example.

Alternative: the `finess.viz` package
-------------------------------------

A second set of plotting tools lives in the Python package
[`$FINESS/python/finess/viz`](../../python/finess/viz).  It reads the
parameter dump `output/parameters.ini.dump` written by each run and is driven
by the generic scripts in `$FINESS/viz`:

    $ python3 $FINESS/viz/plot1d_generic.py          [output/parameters.ini.dump]
    $ python3 $FINESS/viz/plot2d_generic.py          [output/parameters.ini.dump]
    $ python3 $FINESS/viz/plot2d_generic_contour.py  [output/parameters.ini.dump]

These require `$FINESS/python` to be on `$PYTHONPATH` (see
`$FINESS/util/setenv.py`).  Some applications, for example the 2D MHD
`cloud_shock` and `rotated_alfven/convergence` scripts, use this package
directly.
