# Generate parameters ini files for convergence study


# parameters template. 
finess_data_template = '''
; Parameters common to FINESS applications
[finess]
ndims       = 3          ; 1, 2, or 3
nout        = 1          ; number of output times to print results
tfinal      = 4 
initial_dt  = 0.001        ; initial dt
max_dt      = 1.0e10     ; max allowable dt 
max_cfl     = %(max_cfl)f       ; max allowable Courant number
desired_cfl = %(cfl)f    ; desired Courant number
nv          = 500000     ; max number of time steps per call to DogSolve
time_stepping_method = %(ts_method_str)s ; (e.g., Runge-Kutta, Lax-Wendroff, Multiderivative, User-Defined)
space_order = %(s_order)i   ; order of accuracy in space
time_order  = %(t_order)i   ; order of accuracy in time
verbosity   = 1   ; verbosity of output
mcapa       = 0   ; mcapa (capacity function index in aux arrays)
maux        = 0   ; maux (number of aux arrays, maux >= mcapa)
source_term = false   ; source term (1-yes, 0-no)
meqn        = 8   ; number of equations
output_dir  = %(output)s ; location of the output directory
datafmt     = Silo
global_alpha = true
; -------------------------------------------------
;   Cartesian grid data (ignored if Unstructured)
; -------------------------------------------------
[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
my    =  %(my)i  ; number of grid elements in y-direction
mz    =  %(mz)i  ; number of grid elements in z-direction
mbc   =   %(mbc)i      ; number of ghost cells on each boundary
xlow  =   0.0  ; left end point
xhigh =   1.0  ; right end point
ylow  =   0.0  ; lower end point
yhigh =   1.0  ; upper end point
zlow  =   0.0  ; lower end point
zhigh =   1.0  ; upper end point

; -------------------------------------------------
;   WENO values
; -------------------------------------------------
[weno]
weno_version  = JS     ; type of WENO reconstruction (e.g. JS, FD, Z)
epsilon       = 1e-29   ; regulization parameter  ( epsilon > 0.0        )
alpha_scaling = 1.1    ; scaling parameter       ( alpha_scaling >= 1.0 )

[mhd]
gamma = 1.666666666666666667
'''

def generate_parameters_ini_files():
    nruns = 3 
    max_cfl = 0.6
    cfl = 0.5
    ts_method_str = "Lax-Wendroff"
    s_order = 5
    t_order = 3
    mx = 16 
    my = 16
    mz = 16
    mbc = 6
    for nframe in range(nruns):
        parameters_ini_filename = "parameters%02d.ini" % nframe
        output_dir = "output%02d" % nframe
        parameters_dict = {"cfl": cfl,
                            "max_cfl": max_cfl,
                           "ts_method_str": ts_method_str,
                           "s_order": s_order,
                           "t_order": t_order,
                           "output": output_dir,
                           "mx": mx,
                           "my": my,
                           "mz": mz,
                           "mbc": mbc  }
        with open(parameters_ini_filename, 'w') as f:
            f.write(finess_data_template % parameters_dict)
        mx *= 2
        my *= 2
        mz *= 2


if __name__ == "__main__":
    generate_parameters_ini_files()

