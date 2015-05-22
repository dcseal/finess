* DONE Correct setup for grid
    * Check Toth's survey for grid setup
    * Tweak on finess.params.dim2

* TODO Correct boundary values
  For double-periodic boundary condition, there is nothing to change.
    * DONE Try inserting exact values of A3 to boundary
      Result: bad magnetic fields at boundary
    * TODO Try extrapolating A3 based on numerical differential operator and double-periodic boundary conditions on magnetic fields

* Add Constrained Transport
    * TODO Add initial condition for magnetic potential.
    * TODO Runge-Kutta
        * DONE Plumb BeforeStep, ConstructL, UpdateSoln, AfterStep to use
          app redefined version.
        * DONE Plumb to use app-defined AfterQinit
        * STARTED Preliminary version of C.T.
          
    * TODO Lax-Wendroff
