## ------------------------------------------------------------------------- ##
def plotdog1(points_per_dir_in="1",
             outputdir="output",
             qhelpname="qhelp.dat",
             qname="q",
             auxname="a",
             plotq1name="plotq1"):
    """Generic code for plotting FINESS output in matplotlib.

Execute via

    $ python $FINESS/viz/python/plotdog1.py

from an application directory.   For help type

    $ python $FINESS/viz/python/plotdog1.py -h
 
to see a list of options.

    Usage: plotdog1( points_per_dir_in="1", outputdir="output")
 
   points_per_dir = points per direction (spatial dimension)
 
   outputdir = location of output directory
"""


    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    from io_routines import read_params
    from io_routines import read_qfile
    
    TF = os.path.exists(outputdir)
    if TF==False:
        print ""
        print "    Directory not found, outputdir =",outputdir
        print ""
        exit()

    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)
    plotq1_file  = os.path.abspath("plotq1.py")
    local_plotq1 = os.path.exists(plotq1_file)
    if local_plotq1==False:
        from plotq1_default import plotq1
    else:
        from plotq1 import plotq1

    params  = read_params(outputdir)    
    meqn    = int( params['meqn'] )
    maux    = int( params['maux'] )
    nplot   = int( params['nout'] )
    mx      = int( params['mx'  ] )
    xlow    = params['xlow' ]
    xhigh   = params['xhigh']

    # TODO - these are hard-coded here
    meth1  = 1
    mx_old = mx
    dx_old = (xhigh-xlow)/mx

    print("")
    print("       outputdir = %s" % outputdir       )
    print("       sorder    = %i" % meth1           )
    print("")

    # Grid information
    dx = (xhigh-xlow)/float(mx)

    xc = np.zeros(mx,float)
    xc[0] = xlow + 0.5*dx
    for i in range(1,mx):
        xc[i] = xc[i-1] + dx

    xc_old = np.zeros(mx_old,float)
    xc_old[0] = xlow + 0.5*dx_old
    for i in range(1,mx):
        xc_old[i] = xc_old[i-1] + dx_old

    q=-1
    tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
    tmp2 = "".join((tmp1," ) ? "))
    m = raw_input(tmp2)
    print ""

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
            plotq1(m-1, meth1, meqn, mx, time, xc, qsoln, auxsoln)

    plt.ioff()
    print ""
## ------------------------------------------------------------------------- ##
    

#----------------------------------------------------------
if __name__=='__main__':
    """
    If executed at command line prompt, simply run plotfin1
    """

    import optparse
    parser = optparse.OptionParser(
        usage=''' %%prog (-h | [-p POINTS_PER_DIR] [-o OUTPUT_DIRECTORY] )
    
%s''' % plotdog1.__doc__)

    parser.add_option('-p', '--points-per-dir', type='int', default=1, 
                       help='''number of points per cell to be plotted.
                       Default = 1''')
    parser.add_option('-o', '--output-directory', type='string', default='output', 
                       help='''OUTPUT_DIR = output directory where coefficients
                       can be located''')

    opts, args = parser.parse_args()
    plotdog1(opts.points_per_dir, opts.output_directory )

