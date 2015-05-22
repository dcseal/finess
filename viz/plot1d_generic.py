"""
The purpose of this routine is to allow access to a generic (1D) plotting
routine for FINESS.  In order to call this function, make sure that the
environment variable $PYTHONPATH includes $FINESS/python.

To execute this file, from an application directory, simply call

$> python $FINESS/viz/plot1d_generic.py

from the command line.

See also: 
    * plot[2,3]d_generic.py (TODO - write these routines) for other plotting routines.
    * $FINESS/viz/matlab for matlab plotting routines.
"""

from finess.viz.dim1 \
  import ask_which_component_and_which_frame_and_plot_red_dots

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("parameters_ini_filename",
                    nargs = '?',
                    default = 'output/parameters.ini.dump',
                    type = str,
                    help = "default: %(default)s")
args = parser.parse_args()

ask_which_component_and_which_frame_and_plot_red_dots(args.parameters_ini_filename)
