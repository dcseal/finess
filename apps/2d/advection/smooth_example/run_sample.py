#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE

from params_template import dogpack_data_template

def main(cfl, ts_method, space_order, time_order, iterations, mx_start, n_start):

    '''runs time stepping method given by ts_method for 'iterations' number of
    times.  The resolution gets increased by a factor of 'mx_ratio' for each
    iteration.  The starting resolution is given by 'mx_start'.   Output is
    written to folders 'output0000n' for n=0...iterations-1.
'''

    data_file = 'parameters.ini'
    ratio = 2

    integrators   = ['Runge-Kutta', 'SDC', 'Lax-Wendroff', 'Multiderivative', 'User-Defined']
    ts_method_str = integrators[ts_method]
    print(ts_method_str)

    my_start = mx_start
    print space_order
    for i in range(iterations):
        mx_now = mx_start * ratio**i
        my_now = my_start * ratio**i
        #dt_now = dt_start / ratio**i

        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        with closing(open(data_file,'w')) as data:
            # print >> data, dogpack_data_template % locals() 
            ## if we had used same names for everything
            my_dictionary = {'s_order' : space_order, 't_order' : time_order,
                    'cfl' : cfl,
                    'mx' : mx_now, 'my' : my_now,
                    "i_now": (i+n_start), 'ts_method_str' : ts_method_str }
            print >> data, dogpack_data_template % my_dictionary

        # if you want to capture script output, do
        #   Popen(thing to run, shell=True, stdout=PIPE).communicate()[0]
        # cmd = './dog.exe -o outputSL%(s_order)i_%(t_order)i_%(i_now)03i' % my_dictionary
        cmd = './dog.exe -o output_%(i_now)03i' % my_dictionary
        print cmd
        call(cmd, shell=True)
        print ''' 
//////////////////////// finished running output directory output%03i //////////
''' % (i+n_start)

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
        args.order[0], args.order[1], args.frames[1], args.frames[0], 0 )

    """

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h | [-i ITERATIONS] [-x MX_START] [-w MY_START] 
                        [-t DT_START] [-s SPACE_ORDER] [-n N_START])
    
%s''' % main.__doc__)

    parser.add_option('-i', '--iterations', type='int', default=6, 
                       help='''number of times we wish to run computation.
                       Default = 6''')
    parser.add_option('-x', '--mx-start', type='int', default=4, 
                       help='''MX_START = number of grid points to start
                       computation with''')
    parser.add_option('-t', '--dt-start', type='float', default=0.44,
                       help='''DT_START = starting dt value''')
    parser.add_option('-s', '--space-order', type='int', default=5,
                       help='''SPACE_ORDER = spatial order of accuracy''')
    parser.add_option('-p', '--time-order', type='int', default=2,
                       help='''TIME_ORDER = time order of accuracy''')
    parser.add_option('-n', '--n-start', type='int', default=0, 
                       help='''N_START = folder number to start output from.
                       Folders are named output000i with 
                       i = N_START...N_START- 1.  The default value is 0''')
    opts, args = parser.parse_args()

#    check that the user provided enough arguments...
#    if not opts.infile or not opts.outfile:
#        parser.error('Both options -i and -o are required. Try -h for help.')
    main(opts.space_order, opts.time_order, opts.iterations, opts.mx_start,
        opts.dt_start, opts.n_start)
"""
