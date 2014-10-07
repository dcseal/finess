* DONE Correct setup for grid
    * Check Toth's survey for grid setup
    * Tweak on finess.params.dim2

* DONE Correct boundary values
  For double-periodic boundary condition, there is nothing to change.

* Add Constrained Transport
    * TODO Add initial condition for magnetic potential.
    * TODO Runge-Kutta
        * TODO Plumb BeforeStep, ConstructL, UpdateSoln, AfterStep to use
          app redefined version.
    * TODO Lax-Wendroff
