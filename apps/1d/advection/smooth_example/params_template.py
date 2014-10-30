# parameters template. 
#
# This should agree with the current params.ini with the exception of the
# fields we wish to modify.
#
# This template is called by run_sample.py for producing many outputs that are
# needed to run a convergence study.
#

finess_data_template = '''
; Parameters common to FINESS applications
[finess]
ndims       = 1         ; 1, 2, or 3
nout        = 1         ; number of output times to print results
tfinal      = 1.0       ; final time
initial_dt  = 1.0       ; initial dt
max_dt      = 1.0e10    ; max allowable dt 
max_cfl     = 1.00      ; max allowable Courant number
desired_cfl = %(cfl)f   ; desired Courant number
nv          = 500000    ; max number of time steps per call to DogSolve
time_stepping_method = %(ts_method_str)s ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
space_order = %(s_order)i   ; order of accuracy in space
time_order  = %(t_order)i   ; order of accuracy in time
meqn        = 1      ; number of equations
verbosity   = 1      ; verbosity of output
mcapa       = 0      ; mcapa (capacity function index in aux arrays)
maux        = 1      ; maux (number of aux arrays, maux >= mcapa)
source_term = false  ; source term (true or false)
output_dir  = %(output)s ; location of the output directory

[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
mbc   = 5        ; number of ghost cells on each boundary
xlow  =   0.0e0  ; left end point
xhigh =   1.0e0  ; right end point
'''
