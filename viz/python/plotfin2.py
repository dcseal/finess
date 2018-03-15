def plotfin2( outputdir="output", qname="q", auxname="a", plotq2name="plotq2"):
    """ Generic code for plotting FINESS output in matplotlib.

Execute via

    $ python $FINESS/viz/python/plotfin2.py

from an application directory.   For help type

    $ python $FINESS/viz/python/plotfin2.py -h
 
to see a list of options.

    Usage: plotfin2( points_per_dir_in="1", outputdir="output", point_type_in="1", qname="q")
 
    points_per_dir = points per direction (spatial dimension)
 
    outputdir = location of output directory

    point_type = 1:   uniform points on each element
               = 2:   Gauss-Legendre points on each element

    qname          = name of variable filename (i.e., [qname]0000.dat, [qname]0001.dat, ...)
    
"""

    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    from io_routines import read_qfile                  # Reading data files
    from io_routines import parse_ini_parameters        # *.ini parser
   
    TF = os.path.exists(outputdir)
    if(TF==False):
        print("\n    Directory not found, outputdir = %s\n" % outputdir )
        exit()

    # Pull the personalized plotting routine, plotq2.py
    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)
    plotq2_file  = os.path.abspath("plotq2.py")
    local_plotq2 = os.path.exists(plotq2_file)
    if( local_plotq2==False ):
        from plotq2_default import plotq2
    else:
        from plotq2 import plotq2

    ini_params   = parse_ini_parameters( outputdir+'/parameters.ini' )

    # Additional information
    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    nplot    = ini_params['nout']

    # Grid information
    mx      = ini_params['mx'   ]
    xlow    = ini_params['xlow' ]
    xhigh   = ini_params['xhigh']

    my      = ini_params['my'   ]
    ylow    = ini_params['ylow' ]
    yhigh   = ini_params['yhigh']

    # Derived parameters 
    # (TODO - in the future I would like to merge 1D and 2D
    # plotting routines - the outside interface is identical for both! -DS 9/8/2015)
    dx     = (xhigh-xlow)/mx
    dy     = (yhigh-ylow)/my

    print("")
    print("       outputdir = %s" % outputdir       )
    print("")

    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)

    # Set up the grid
    xc = np.zeros((mx,my),float)
    yc = np.zeros((mx,my),float)
    for j in range(0,my+1):
        xl[:,j] = xlow + dx*np.arange(mx+1)[:]
    for i in range(0,mx+1):
        yl[i,:] = ylow + dy*np.arange(my+1)[:]

    for j in range(0,my):
        xc[:,j] = (xlow+0.5*dx) + dx*np.arange(mx)[:]
    for i in range(0,mx):
        yc[i,:] = (ylow+0.5*dy) + dy*np.arange(my)[:]

    q=-1;
    tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
    tmp2 = "".join((tmp1," ) ? "))
    m = raw_input(tmp2)
    print("")

    # Select an equation to plot. Defuaul: m=1
    if (not m):
        m = 1
    else:
        m = int(m)

    if m<1:
        print ""
        print "  Error, need m > 1,  m = ",m
        print ""
        return -1

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

        if (nf!=-1):

            ## -------------------------------------------------------------- ##
            # Solution -- q
            # solution should be found in file
            #     outputdir/q[n1].dat
            ## -------------------------------------------------------------- ##
            qfile_tmp_tmp = "".join((str(n1+10000),".dat"))
            qfile_tmp = qname + qfile_tmp_tmp[1:]
            qfile = "".join(("".join((outputdir,"/")),qfile_tmp))

            qtmp    = np.zeros(mx*my*meqn,float)   
            time    = read_qfile(mx*my*meqn,qfile,qtmp)
            qcoeffs = np.reshape(qtmp,(meqn,my,mx))
            
            qsoln = np.zeros((mx*my,meqn),float)
            sample_state2_cart_mod(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln)
            qaug = np.zeros((mx+1,my+1,meqn),float)
            qaug[0:mx,0:my,0:meqn] = np.reshape(qsoln,(mx,my,meqn),'F')

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
            plotq2(xc, yc, qsoln, auxsoln, time, ini_params, m-1 )

    plt.ioff()
    print ""

## ------------------------------------------------------------------------- ##
def parse_input( help_message ):
    """Parse command line arguments for 2D plotting routines."""

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

    parser.add_argument( '-a', '--plotq2-name', 
        type = str, 
        default     = 'plotq2', 
        help        =
'''filename for file with additional plotting options.
(default: plotq2)''')

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
    args = parse_input( plotfin2.__doc__ )

    # Call the main 1D plotting routine
    plotfin2(args.outputdir, args.q_name, args.aux_name, args.plotq2_name )

