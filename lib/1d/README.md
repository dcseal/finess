FINESS - the FINite difference ESSentially non-oscillatory software package.
============================================================================

This folder contains all of the main library routines required for running
[FINESS](../..//README.md).

Basic structure
---------------

The most expensive part of the code is the function

    ConstructL

which defines the right hand side of the Method of Lines (MOL) formulation of
the PDE: 

    q_t = L( q ).

Required functions and templates
--------------------------------

See 

    $FINESS/lib/1d/blanks

for required functions that need to be implemented by the user in order to run
the code.  These include the following routines:

* FluxFunc          - flux function for the PDE.
* QinitFunc         - initial conditions for the PDE.
* SourceTermFunc    - source term function, if any.  (set ... = 0 if not used)
* AuxFunc           - auxiliary function (optional, set maux = 0 if not used)
* SetWaveSpd        - compute the largest and smallest wave speed of the PDE.
* SetBndValues      - boundary value data
* ProjectLeftEig    - Projection onto the left eigenvalues of the system
* ProjectRightEig   - ... right eigenvalues of the system

It is recommended that a user copies an existing application and then modifies
it in order to implement a new solver.

WENO reconstruction options
===========================

TODO - currently no WENO reconstruction options are implemented.
