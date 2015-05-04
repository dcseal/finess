# Introduction #

[FINESS](./README.md) stands for the FINite difference ESSentially non-oscillatory software
package.  It is a set of C++ libraries intended to be used by applications to
solve hyperbolic partial differential equations.  This code was heavily
derived from DoGPack (http://www.dogpack-code.org).

See: [$FINESS/LICENSE.txt](LICENSE.txt) for the end user license agreement.

## Overview ##

A typical FINESS application solves a hyperbolic equation of the form 

       q_t + div(f) = \psi.

The user is expected to supply the flux function f and the source term \psi as
library callbacks, as defined in the [structure of code](#Structure_of_code)
section of this document.  Plotting options can be found in the
[plotting section](#plotting) of this document.

The installation procedure is essentially idential to that of 
[DoGPack](http://www.dogpack-code.org/install.html), where you essentially
replace each occurance of $DOGPACK with $FINESS.  The public source code can
be pull by calling

    $ git clone https://bitbucket.org/dseal/finess-release

from the command line.

<a name="Structure_of_code"></a>
## Structure of code ##

Every application lives in 

    $FINESS/apps/

This directory contains problem definitions including initial conditions and
flux functions.  The main library is located in 

    $FINESS/lib/

which defines various time stepping routines, WENO reconstructions and the
main control flow structures.  Further dimension specific libraries are
located in

    $FINESS/lib/1d/
    $FINESS/lib/2d/

Some documentation for the [1d](lib/1d/README.md) and [2d](lib/2d/README.md)
has been written, but most of it lives inside the source code itself.
Most applications interface with the library by relinking the source files in
their relevant Makefiles.  More on this later.

## How to compile and run an example ##

The best way to learn how to run the code is to compile and run a working
example.  This section defines how to compile and run a single 1D example,
which is located here:

    $FINESS/apps/1d/advection/smooth_example

But first, what's up with the "$FINESS" environment variable?  Well, you need
to set that up.

### Setting up the environment ###

NOTE: these steps are identical for setting up the environment variable for
[DoGPack](http://www.dogpack-code.org/install.html).

The easiest way to set up the $FINESS environment variable is to use the short
script provided in the util subdirectory.  Go into your base directory you
created above and run

    $ python util/setenv.py

This script should produce two files that contain the shell script for setting
the above variables. By default these files are called setenv.bash and
setenv.csh. These can be used by running

    $ source setenv.bash

or

    $ source setenv.csh

depending on your shell (this can be checked by typing 'echo $SHELL') from
your command line). It is highly recommended that you place these in your 
.bashrc, .cshrc, or .profile file to be run automatically when you open a
terminal, otherwise you will have to the step of 'sourcing' these files every
time you open a new terminal.

* If using MATLAB for plotting, you may need to manually set the path in MATLAB
using the MATLAB command window. The correct path is shown in the setenv.bash
and setenv.csh files.

### Testing the installation ###

Once the installation is completed, the next step is to make sure that the
installation and setup were successful. Do this by entering a particular
example and compiling the code by typing:

    $ cd $FINESS/apps/1d/advection/smooth_example
    $ make

Once compiled, execute the code by typing:

    $ finess.exe

You should see a bunch output printed to the screen, followed by a new folder
called 'output' that the code should have created to store data for the
simulation.  To redirect output to a different folder, you can replace the
"output_folder" argument inside the parameters.ini file.  For example, if you
insert the line

    output_dir  = output ; location of the output directory

with 

    output_dir  = some_other_folder ; location of the output directory

then the code will create

    $FINESS/apps/1d/advection/smooth_example/some_other_folder

and save data there.  Note that the parameters.ini file used to create this
run is saved, so you can reproduce results from this run by copying that file
back to the application directory.

<a name="Plotting"></a>
## Plotting options ##

After you run an application, you will likely want to plot the results.
[FINESS](README.md) supports both [MATLAB](viz/matlab/README.md) and
[Python](viz/python/README.md) plotting routines.  All of these can be found in
the

    $FINESS/viz/

directory.  The basic format requires the user to call the main plotting
routines from an application directory.  Local options, such as grid labels,
axis options etc. are supplied by a single script defined by each application.


