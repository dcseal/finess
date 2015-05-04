# parameters template. 
#
# This should agree with the current params.ini with the exception of the
# fields we wish to modify.
#
# This template is called by run_sample.py for producing many outputs that are
# needed to run a convergence study.
#

dogpack_data_template = '''
; Parameters common to dogpack applications
[dogParams]
defaults_file = "$FINESS/config/dogParams_defaults.ini"
ndims       = 2          ; 1 or 2
mesh_type   = Cartesian  ; (either Cartesian or Unstructured) 
nout        = 1          ; number of output times to print results
tfinal      = 2.0        ; final time
dtv(1)      = 1.0e0     ; initial dt
dtv(2)      = 1.1e0      ; max allowable dt 
cflv(1)     = 1.00       ; max allowable Courant number
cflv(2)     = %(cfl)f       ; desired Courant number
nv          = 500000     ; max number of time steps per call to DogSolve
time_stepping_method = %(ts_method_str)s ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
limiter_method = moment ; (e.g., moment, viscosity)
space_order = %(s_order)i   ; =method(1)= order of accuracy in space
time_order  = %(t_order)i   ; =method(2)= order of accuracy in time
use_limiter = 0   ; =method(3)= use limiter (1-yes, 0-no)
verbosity   = 1   ; =method(4)= verbosity of output
mcapa       = 0   ; =method(5)= mcapa (capacity function index in aux arrays)
maux        = 0   ; =method(6)= maux (number of aux arrays, maux >= mcapa)
source_term = 0   ; =method(7)= source term (1-yes, 0-no)
meqn        = 5   ; number of equations
mrestart    = 0   ; restart from old data (1-yes, 0-no)
nstart      = 0   ; if mrestart==1: from file q(nstart)_restart.dat
datafmt     = 1   ; 1 for ascii, 5 for hdf5.
ic_quad_order = 20  ; order of quadrature for L2-projection of initial conditions
; -------------------------------------------------
;   Cartesian grid data (ignored if Unstructured)
; -------------------------------------------------
[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
my    =  %(my)i  ; number of grid elements in y-direction
mbc   =   5      ; number of ghost cells on each boundary
xlow  =   0.0e0  ; left end point
xhigh =   2.0e0  ; right end point
ylow  =   0.0e0  ; lower end point
yhigh =   2.0e0  ; upper end point

[weno]
weno_version  = JS     ; type of WENO reconstruction (e.g. JS, FD, Z)
epsilon       = 1e-29  ; regulization parameter  ( epsilon > 0.0        )
alpha_scaling = 1.0    ; scaling parameter       ( alpha_scaling >= 1.0 )
'''


