
# parameters template. 
#
# This should agree with the current parameters.ini with the exception of the
# fields we wish to modify.
finess_data_template = '''
; Parameters common to FINESS applications
[finess]
ndims       = 2          ; 1, 2, or 3
nout        = 1          ; number of output times to print results
;tfinal      = 2.0        ; final time
tfinal      = 0.80      ; final time (for RP1)
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
maux        = 0   ; maux (number of aux arrays, maux >= mcapa)
source_term = false   ; source term (1-yes, 0-no)
meqn        = 5   ; number of equations
output_dir  = %(output)s ; location of the output directory

; -------------------------------------------------
;   Cartesian grid data (ignored if Unstructured)
; -------------------------------------------------
[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
my    =  %(my)i  ; number of grid elements in y-direction
mbc   =   11      ; number of ghost cells on each boundary

;xlow  =  0.0e0  ; left end point
;xhigh =  1.0e0  ; right end point
;ylow  =  0.0e0  ; lower end point
;yhigh =  1.0e0  ; upper end point

; Grid spacing for RP1 and RP3
xlow  = -0.5e0  ; left end point
xhigh =  0.5e0  ; right end point
ylow  = -0.5e0  ; lower end point
yhigh =  0.5e0  ; upper end point

; -------------------------------------------------
;   WENO values
; -------------------------------------------------
[weno]
weno_version  = JS     ; type of WENO reconstruction (e.g. JS, FD, Z)
epsilon       = 1e-29  ; regulization parameter  ( epsilon > 0.0        )
alpha_scaling = 1.1    ; scaling parameter       ( alpha_scaling >= 1.0 )

; -------------------------------------------------
;  Euler parameters
; -------------------------------------------------
[euler]
gamma = 1.4                 ; gas constant
riemann_problem_number = 1  ; Option [0,1,2,3,4,5] for which Riemann problem to solver.  
'''

qsub_template = '''
#!/bin/bash -login
#PBS -l walltime=48:02:00,nodes=1:ppn=%(proc_per_node)i
#PBS -j oe
#PBS -N run_2d_euler_%(descriptor)s
 
# load necessary modules, e.g.
# module load HDF5
module load gcc
module load Python
 
# change to the working directory where your code is located
#
# TODO - I don't know how to get this to use the environment variable $(FINESS)
# that should be loaded with each login shell.  I tried doing the following:
#
# source ~/.bashrc
#
# but that didn't seem to work
#
# Instead, I am currently hard coding my own folder to this folder change:
cd ~/code/finess/apps/2d/euler/two_d_riemann_test
#cd $(PBS_O_WORKDIR)

# Compile your code (uncomment this if you do not do this by hand)
# make cleanallo; make -j;

export OMP_NUM_THREADS='24'

# call your executable
pwd
aprun -n1 -d24 ./finess.exe %(params_file_name)s 
'''
