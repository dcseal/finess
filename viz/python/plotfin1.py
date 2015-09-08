## ------------------------------------------------------------------------- ##
def plotfin1( outputdir="output", qname="q", auxname="a", plotq1name="plotq1"):
    """Generic code for plotting FINESS output in matplotlib.

Execute via

    $ python $FINESS/viz/python/plotdog1.py

from an application directory.   For help type

    $ python $FINESS/viz/python/plotdog1.py -h
 
to see a list of options.

    Usage: plotdog1( outputdir="output")
 
   points_per_dir = points per direction (spatial dimension)
 
   outputdir = location of output directory
"""


    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    from io_routines import parse_ini_parameters            
    from io_routines import read_qfile                  # Reading data files
    from io_routines import parse_ini_parameters        # *.ini parser

    TF = os.path.exists( outputdir )
    if TF==False:
        print("\n    Directory not found, outputdir = %s\n" % outputdir )
        exit()

    # Pull the personalized plotting routine, plotq1.py
    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)
    plotq1_file  = os.path.abspath("plotq1.py")
    local_plotq1 = os.path.exists(plotq1_file)
    if( local_plotq1==False ):
        from plotq1_default import plotq1
    else:
        from plotq1 import plotq1

    ini_params   = parse_ini_parameters( outputdir+'/parameters.ini' )

    mx      = ini_params['mx'   ]
    xlow    = ini_params['xlow' ]
    xhigh   = ini_params['xhigh']

    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    nplot    = ini_params['nout']
    mx       = ini_params['mx']
    xlow     = ini_params['xlow']
    xhigh    = ini_params['xhigh']

    # TODO - these are hard-coded here
    # Derived parameters:
    dx     = (xhigh-xlow)/mx

    print("")
    print("       outputdir = %s" % outputdir       )
    print("")

    # Grid information
    dx = (xhigh-xlow)/float(mx)
    xc = np.linspace( xlow+0.5*dx, xhigh-0.5*dx, mx )

    q=-1
    tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
    tmp2 = "".join((tmp1," ) ? "))
    m = raw_input(tmp2)
    print("")

    # Select an equation to plot. Defuaul: m=1
    if(not m):
        m = 1
    else:
        m = int(m)

    # Check user input
    if( m < 1 or m > meqn ):
        print("\n  Error, need 1 <= m <= meqn = %d.  You chose m = %d \n" % (meqn, m) )
        exit(1)

    n  = 0;
    nf = 0;
    n1 = -1;

    plt.ion()
    while (nf!=-1):

        tmp1 = "".join((" Plot which frame ( 0 - ",str(nplot)))
        tmp2 = "".join((tmp1," ) [type -1 or q to quit] ? "))
        nf = raw_input(tmp2)
        if (not nf):
            n1 = n1 + 1
            nf = 0
        elif(nf=="q"):
            nf = -1
        else:
            nf = int(nf)
            n1 = nf

        if(n1>nplot):
            print("\n End of plots \n")
            n1 = nplot

        if(nf!=-1):

            ## -------------------------------------------------------------- ##
            # Solution -- q
            # solution should be found in file
            #     outputdir/q[n1].dat
            ## -------------------------------------------------------------- ##
            qfile_tmp_tmp = "".join((str(n1+10000),".dat"))
            qfile_tmp     = "q" + qfile_tmp_tmp[1:]
            qfile         = "".join(("".join((outputdir,"/")),qfile_tmp))

            time,qsoln  = read_qfile(mx, meqn, qfile )

            if(maux>0):

                ## ---------------------------------------------------------- ##
                # Aux arrays -- aux
                # solution should be found in file
                #     outputdir/aux[n1].dat
                ## ---------------------------------------------------------- ##
                afile_tmp_tmp = "".join((str(n1+10000),".dat"))
                afile_tmp = auxname + afile_tmp_tmp[1:]
                afile = "".join(("".join((outputdir,"/")),afile_tmp))

                time_aux, auxsoln = read_qfile(mx, maux, afile )

            else:
                auxsoln = 0.0;

            # USER SUPPLIED FUNCTION (or default function )
            plotq1(xc, qsoln, auxsoln, time, ini_params, m-1 )

    plt.ioff()
    print ""
## ------------------------------------------------------------------------- ##
    
## ------------------------------------------------------------------------- ##
def parse_input( help_message ):
    """Parse command line arguments for 1D plotting routines."""

    import argparse, sys

    parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description =    help_message,
      formatter_class = argparse.RawTextHelpFormatter,
    )

    parser.add_argument('-o', '--outputdir', 
        type        = str, 
        default     = 'output', 
        help        =
'''Location of the output directory where coefficients can be located.
(default: output)''')

    parser.add_argument( '-q', '--q-name', type = str, default='q', 
        help        =
'''Name of variable used in the filename.  In most routines this is 'q' but in some cases
one may wish to save additional data.  For example, a 2D code may want to save
1D data objects, and then plot them.
(default: q).''')

    parser.add_argument( '-x', '--aux-name', type = str, default='a', 
        help        =
'''Name of aux variable used in the filename.  In most routines this is 'a' but in some cases
one may wish to save additional data.  For example, a 2D code may want to save
1D data objects, and then plot them.
(default: a).''')

    parser.add_argument( '-a', '--plotq1-name', 
        type = str, 
        default     = 'plotq1', 
        help        =
'''filename for file with additional plotting options.
(default: plotq1)''')

    return parser.parse_args()
#----------------------------------------------------------


#----------------------------------------------------------
if __name__== '__main__':
    """Default python plotting routine for 1D simulations in FINESS.

    When run from the command line, this script parses user supplied command
    line arguments, if any, and then executes plotfin1.  To see a list of
    options, type 

        python $FINESS/viz/python/plotfin1.py -h

    """

    # Parse input arguments
    args = parse_input( plotfin1.__doc__ )

    # Call the main 1D plotting routine
    plotfin1(args.outputdir, args.q_name, args.aux_name, args.plotq1_name )

