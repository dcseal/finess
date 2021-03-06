finess_data_template = '''
; Parameters common to FINESS applications
[finess]
ndims       = 1         ; 1, 2, or 3
nout        = 1             ; number of output times to print results
tfinal      = 1.591549430918953e-01 ; final time 
initial_dt  = 1.00e-03   ; initial dt
max_dt      = 1.0e10     ; max allowable dt 
max_cfl     = 0.84       ; max allowable Courant number
desired_cfl = %(cfl)g    ; desired Courant number
nv          = 500000    ; max number of time steps per call to DogSolve
time_stepping_method = %(ts_method_str)s ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
space_order = %(s_order)i   ; =method(1)= order of accuracy in space
time_order  = %(t_order)i   ; =method(2)= order of accuracy in time
verbosity   = 1   ; verbosity of output
mcapa       = 0   ; mcapa (capacity function index in aux arrays)
maux        = 0   ; maux (number of aux arrays, maux >= mcapa)
source_term = false   ; source term (true or false)
meqn        = 1   ; number of equations
output_dir  = %(output)s ; location of the output directory
[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
mbc   =      5   ; number of ghost cells on each boundary
xlow  =  0.0e0   ; left end point
xhigh =  2.0e0   ; right end point
'''
