# FINESS - the FINite difference ESSentially non-oscillatory software package.  #

This folder contains all of the top-level library routines that required for
running of [FINESS](../..//README.md) in any dimension.

For further details on the 1- and 2D library routines, please visit the
corresponding directories (e.g. [$FINESS/lib/1d/](1d/README.md) and
[$FINESS/lib/2d/](2d/README.md) ).

## WENO reconstruction options ##

The WENO reconstruction options are concisely listed in 

* WenoReconstruct.h

A user can specify this option in the [weno] section of any parameters.ini
file, located in each application directory.  For example, a typical section
may look like

    weno_version  = JS     ; type of WENO reconstruction (e.g. JS, FD, Z)
    epsilon       = 1e-29  ; regulization parameter  ( epsilon > 0.0        )
    alpha_scaling = 1.1    ; scaling parameter       ( alpha_scaling >= 1.0 )

The description of these parameters in the source code is located in

    $FINESS/lib/WenoParams.h.

We currently have orders 5, 7, 9 and 11 implented with the Jiang and Shu method,
and orders 5, 7 and 9 implemented for the WENO-Z option.  The WENO
reconstruction is typically called from ConstructL, which depends on the
dimension of the problem.
