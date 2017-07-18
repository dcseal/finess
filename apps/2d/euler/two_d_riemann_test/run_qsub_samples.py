#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE

from params_template import finess_data_template, qsub_template

from time import sleep

def main(ts_method):
    '''Create multiple qsub calls and run the code.'''

    integrators   = ['Runge-Kutta', 'SDC', 'Lax-Wendroff', 'Multiderivative', 'User-Defined']
    ts_method_str = integrators[ts_method]
    print(ts_method_str)

    # Parameters to run samples over.  This script looks at all combinations of
    # these parameters
    cfl_range = [0.3,0.4]
    sorders   = [5,7,9,11]
    mx_range  = [64,128,256,512,728,999]
    mx_range  = [512,728,999]

    cfl_range = [0.3]
    sorders   = [7,9]
    mx_range  = [64,128,256,512,728,999]
    riemann_problems = [2]

    for mx_now in mx_range:
      for space_order in sorders:
        for riemann_prob_num in riemann_problems:
            for cfl in cfl_range:

                mx_now = mx_now
                my_now = mx_now

                if( ts_method == 0 ):
                    time_order = 4
                else:
                    time_order = space_order

                my_dictionary = {'s_order' : space_order, 't_order' : time_order,
                        'cfl' : cfl, 'mx' : mx_now, 'my' : my_now,
                        'rp_num' : riemann_prob_num,
                        'ts_method_str' : ts_method_str }

                if( ts_method == 0 ):
                    my_dictionary['output'] = 'output%(mx)03i_rk_%(s_order)02i_cfl%(cfl)2.1f' % my_dictionary
                    params_file_name = 'parameters%(mx)03i_rk_%(s_order)02i_cfl_%(cfl)2.1f.ini' % my_dictionary
                else:
                    my_dictionary['output'] = 'output%(mx)03i_dt_%(s_order)02i_cfl%(cfl)2.1f' % my_dictionary
                    params_file_name = 'parameters%(mx)03i_dt_%(s_order)02i_cfl_%(cfl)2.1f.ini' % my_dictionary

                with closing(open(params_file_name,'w')) as data:
                    print >> data, finess_data_template% my_dictionary


                my_dictionary['proc_per_node'] = 24


                if( ts_method == 0 ):
                    my_dictionary['descriptor'] = '%(mx)03i_rk_%(s_order)02i_cfl_%(cfl)2.1f' % my_dictionary
                    my_dictionary['params_file_name'] = params_file_name
                    qsub_file_name = 'run_app_%(mx)03i_rk_%(s_order)02i_cfl_%(cfl)2.1f.qsub' % my_dictionary
                else:
                    my_dictionary['descriptor'] = '%(mx)03i_dt_%(s_order)02i_cfl_%(cfl)2.1f' % my_dictionary
                    my_dictionary['params_file_name'] = params_file_name
                    qsub_file_name = 'run_app_%(mx)03i_dt_%(s_order)02i_cfl_%(cfl)2.1f.qsub' % my_dictionary

                with closing(open(qsub_file_name,'w')) as data:
                    print >> data, qsub_template % my_dictionary

                # if you want to capture script output, do
                #   Popen(thing to run, shell=True, stdout=PIPE).communicate()[0]
                cmd = 'qsub %s' % qsub_file_name
                print('Running command', cmd)
                call(cmd, shell=True)
                sleep(1)

def parse_input( help_message ):
    """Parse input arguments."""
  
    import argparse, sys
  
    parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description =    help_message,
      formatter_class = argparse.RawTextHelpFormatter,
    )

    """
    parser.add_argument('CFL',
                      type = float,
                      help = 'maximum Courant number in domain')
 
    """

    parser.add_argument('-s','--time_integrator',
                      type    = int,
                      choices = range(5),
                      default =  4,
                      dest    = 't_stepper',
                      metavar = 'X',
                      help    = 
    '''choose integrator X for time-stepping:
  0. Runge-Kutta
  1. SDC
  2. Lax-Wendroff
  3. Multiderivative
  4. 'User-Defined' time integrator.  (See Makefile for what gets linked to)
(default: 4)''')

    return parser.parse_args()

if __name__ == '__main__':
    """ Batch refinement study. 
    
    This script runs multiple samples, which can be later analyzed for a
    convergence study."""

    # Parse input arguments
    args = parse_input( main.__doc__ )
    print(args)
    print('')

    main( args.t_stepper )
