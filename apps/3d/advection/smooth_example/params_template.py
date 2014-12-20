# parameters template. 
#
# This should agree with the current params.ini with the exception of the
# fields we wish to modify.
#
# This template is called by run_sample.py for producing many outputs that are
# needed to run a convergence study.

finess_data_template = '''
; Parameters common to FINESS applications
[finess]
ndims       = 3          ; 1, 2, or 3
nout        = 1          ; number of output times to print results
tfinal      = 1.0        ; final time
initial_dt  = 1.0        ; initial dt
max_dt      = 1.0e10     ; max allowable dt 
max_cfl     = 1.00       ; max allowable Courant number
desired_cfl = %(cfl)f    ; desired Courant number
nv          = 500000     ; max number of time steps per call to DogSolve
time_stepping_method = %(ts_method_str)s ; (e.g., Runge-Kutta, Lax-Wendroff, Multiderivative, User-Defined)
space_order = %(s_order)i   ; order of accuracy in space
time_order  = %(t_order)i   ; order of accuracy in time
verbosity   = 1   ; verbosity of output
mcapa       = 0   ; mcapa (capacity function index in aux arrays)
maux        = 3   ; maux (number of aux arrays, maux >= mcapa)
source_term = false   ; source term (1-yes, 0-no)
meqn        = 1   ; number of equations
output_dir  = %(output)s ; location of the output directory

; -------------------------------------------------
;   Cartesian grid data (ignored if Unstructured)
; -------------------------------------------------
[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
my    =  %(my)i  ; number of grid elements in y-direction
mz    =  %(mz)i  ; number of grid elements in z-direction
mbc   =   3      ; number of ghost cells on each boundary
xlow  =   0.0e0  ; left end point
xhigh =   1.0e0  ; right end point
ylow  =   0.0e0  ; lower end point
yhigh =   1.0e0  ; upper end point
zlow  =   0.0e0  ; lower end point
zhigh =   1.0e0  ; upper end point

; -------------------------------------------------
;   WENO values
; -------------------------------------------------
[weno]
weno_version  = FD     ; type of WENO reconstruction (e.g. JS, FD, Z)
epsilon       = 1e-6   ; regulization parameter  ( epsilon > 0.0        )
alpha_scaling = 1.0    ; scaling parameter       ( alpha_scaling >= 1.0 )
'''
