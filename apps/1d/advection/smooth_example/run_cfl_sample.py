#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE

dogpack_data_template = '''
; Parameters common to dogpack applications
[dogParams]
defaults_file = "$FINESS/config/dogParams_defaults.ini"
ndims       = 1             ; 1 or 2
mesh_type   = Cartesian     ; (either Cartesian or Unstructured) 
nout        = 1             ; number of output times to print results
tfinal      = %(t_final)f   ; final time
dtv(1)      = 1.0           ; initial dt
dtv(2)      = 1e10          ; max allowable dt 
cflv(1)     = 1.0           ; max allowable Courant number
cflv(2)     = %(cfl)f       ; desired Courant number
nv          = 500000        ; max number of time steps per call to DogSolve
time_stepping_method = %(ts_method_str)s ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
limiter_method = moment ; (e.g., moment, viscosity)
space_order = %(s_order)i   ; =method(1)= order of accuracy in space
time_order  = %(t_order)i   ; =method(2)= order of accuracy in time
use_limiter = 0   ; =method(3)= use limiter (1-yes, 0-no)
verbosity   = 0   ; =method(4)= verbosity of output
mcapa       = 0   ; =method(5)= mcapa (capacity function index in aux arrays)
maux        = 1   ; =method(6)= maux (number of aux arrays, maux >= mcapa)
source_term = 0   ; =method(7)= source term (1-yes, 0-no)
meqn        = 1   ; number of equations
mrestart    = 0   ; restart from old data (1-yes, 0-no)
nstart      = 0   ; if mrestart==1: from file q(nstart)_restart.dat
datafmt     = 1   ; 1 for ascii, 5 for hdf5.
withPyClawPlotting = 0; (1-yes, 0-no)

[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
mbc   =      5   ; number of ghost cells on each boundary
xlow  = -1.0e0   ; left end point
xhigh =  1.0e0   ; right end point

[weno]
weno_version  = FD     ; type of WENO reconstruction (e.g. JS, FD, Z)
;weno_version  = JS     ; type of WENO reconstruction (e.g. JS, FD, Z)
;weno_version  = Z      ; type of WENO reconstruction (e.g. JS, FD, Z)
epsilon       = 1e-29  ; regulization parameter  ( epsilon > 0.0        )
alpha_scaling = 1.0    ; scaling parameter       ( alpha_scaling >= 1.0 )
'''
    
def main(cfl_vec, num_frames, ts_method, space_order, time_order, mx_start, n_start, t_final ):
    '''Simple script for performing a CFL scan in order to compute total
    variation of the problem.
'''
    data_file = 'parameters.ini'
    ratio = 2

    integrators   = ['Runge-Kutta', 'SDC', 'Lax-Wendroff', 'Multiderivative', 'User-Defined']
    ts_method_str = integrators[ts_method]
    print(ts_method_str)
    print(cfl_vec)
    print(num_frames)

    dc = (cfl_vec[1] - cfl_vec[0])/(num_frames)

    for i in range( num_frames+1 ):

        with closing(open(data_file,'w')) as data:

            # print >> data, dogpack_data_template % locals() 
            ## if we had used same names for everything

            my_dictionary = {'s_order' : space_order, 't_order' : time_order, 
                    'cfl' : cfl_vec[0] + i*dc,
                    'mx' : mx_start,
                    "i_now": (i+n_start),
                    'ts_method_str' : ts_method_str,
                    't_final' : t_final }
            print >> data, dogpack_data_template % my_dictionary

        # if you want to capture script output, do
        #   Popen(thing to run, shell=True, stdout=PIPE).communicate()[0]
        #%cmd = './dog.exe -o outputSL%(s_order)i_%(t_order)i_%(i_now)03i' % my_dictionary
        cmd = './dog.exe -o output%04i' % my_dictionary['i_now']
        print(cmd)
        call(cmd, shell=True)
        print(''' 
//////////////////////// finished running output directory output%04i //////////
''' % my_dictionary['i_now'] )

def parse_input( help_message ):
  
    import argparse, sys
  
    parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description =    help_message,
      formatter_class = argparse.RawTextHelpFormatter,
    )

    parser.add_argument('CFL_RANGE',
                      type  = float,
                      nargs = 2,
                      help  = '''Range of CFL numbers to search''')

    parser.add_argument('-n','--num_frames',
                      type    = int,
                      default = 100,
                      metavar = ('NFRAMES'),
                      help    = 
''' Number of frames to print.  See also: cfl_range.
(default: 100 )''')

    parser.add_argument('-s','--time_integrator',
                      type    = int,
                      choices = range(4),
                      default =  0,
                      dest    = 't_stepper',
                      metavar = 'X',
                      help    = 
    '''choose integrator X for time-stepping:
  0. Runge-Kutta
  1. SDC
  2. Lax-Wendroff
  3. Multiderivative
  4. 'User-Defined' time integrator.  (See Makefile for what gets linked to)
(default: 0)''')

    parser.add_argument('-m','--mx',
                      type    = int,
                      default = 200,
                      metavar = 'MX',
                      help    = 
''' Number of grid points.
(default: 200 )''')

    parser.add_argument('-o', '--order',
                      type = int,
                      nargs   = 2,
                      default = [5, 3],
                      metavar = ('S_ORDER', 'T_ORDER'),
                      help = 
'''Order of accuracy in space S_ORDER, and time T_ORDER.
(default: [5,3])''')
 
    parser.add_argument('-t', '--t_final',
                        type = float,
                        default = 0.2,
                        metavar = 'T_FINAL',
                        help =
''' Final time of the simulation.
(default: 0.2 )''')

    return parser.parse_args()

if __name__ == '__main__':
    """ Batch refinement study. 
    
    This script runs multiple samples, which can be later analyzed for a
    convergence study."""

    # Parse input arguments
    args = parse_input( main.__doc__ )
    print(args)
    print('')

    main( args.CFL_RANGE, args.num_frames, args.t_stepper, args.order[0], args.order[1], args.mx, 0, args.t_final )
