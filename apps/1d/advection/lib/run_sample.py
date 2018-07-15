#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE

def main(cfl, ts_method, space_order, time_order, iterations, mx_start, n_start, ratio):
    '''Simple script for performing a batch refinement study.

  This script runs time stepping method given by ts_method for 'iterations' number of
  times.  The resolution gets increased by a factor of 'mx_ratio' for each
  iteration.  The starting resolution is given by 'mx_start'.   Output is
  written to folders 'output0000n' for n=0...iterations-1.
'''

    # Need these calls because this script is not in the same directory as the
    # finess_data_template
    import sys,os
    sys.path.append(os.getcwd())
    from params_template import finess_data_template

    data_file = 'parameters.ini'

    integrators   = ['Runge-Kutta', 'SDC', 'Lax-Wendroff', 'Multiderivative', 'User-Defined']
    ts_method_str = integrators[ts_method]
    print(ts_method_str)

    mx_now = mx_start
    for i in range(iterations):

        with closing(open(data_file,'w')) as data:

            ## if we had used same names for everything
            my_dictionary = {'s_order' : space_order, 't_order' : time_order,
                    'cfl' : cfl,
                    'mx' : mx_now, 
                    "i_now": (i+n_start), 'ts_method_str' : ts_method_str,
            }
            my_dictionary['output'] = 'output_%(i_now)04i' % my_dictionary
            print >> data, finess_data_template % my_dictionary

        # if you want to capture script output, do
        #   Popen(thing to run, shell=True, stdout=PIPE).communicate()[0]
        #%cmd = './dog.exe -o outputSL%(s_order)i_%(t_order)i_%(i_now)03i' % my_dictionary
        cmd = './finess.exe'
        print(cmd)
        call(cmd, shell=True)
        print(''' 
//////////////////////// finished running output directory output%03i //////////
''' % (i+n_start) )

        mx_now = int( mx_now*ratio )

def parse_input( help_message ):
  
    import argparse, sys
  
    parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description =    help_message,
      formatter_class = argparse.RawTextHelpFormatter,
    )

    parser.add_argument('CFL',
                      type = float,
                      help = 'maximum Courant number in domain')
  
    parser.add_argument('-s','--time_integrator',
                      type    = int,
                      choices = range(5),
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

    parser.add_argument('-f','--frames',
                      type    = int,
                      nargs   = 2,
                      default = [50, 7],
                      metavar = ('MX_START', 'N_FRAMES'),
                      help    = 
''' Refinment parameters:
Produce N_FRAMES of refinements, starting with 
MX_START grid points.
(default: [50, 7] )''')

    parser.add_argument('-r','--rratio',
                      type    = float,
                      default = 2.0,
                      metavar = 'R_RATIO',
                      help    = 
''' Refinment ratio:
Each refinement for the solution increases the ratio by a factor of R_RATIO.
(default: 2''')


    parser.add_argument('-t', '--order',
                      type = int,
                      nargs   = 2,
                      default = [5, 3],
                      metavar = ('S_ORDER', 'T_ORDER'),
                      help = 
'''Order of accuracy in space S_ORDER, and time T_ORDER.
(default: [5,3])''')
  
    return parser.parse_args()

if __name__ == '__main__':
    """ Batch refinement study. 
    
    This script runs multiple samples, which can be later analyzed for a
    convergence study."""

    # Parse input arguments
    args = parse_input( main.__doc__ )
    print(args)
    print('')

    main( args.CFL, args.t_stepper, 
        args.order[0], args.order[1], args.frames[1], args.frames[0], 0, args.rratio )
