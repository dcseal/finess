## ------------------------------------------------------------------------- ##
def plotdog1(points_per_dir_in="1",
             outputdir="output",
             point_type=1,
             qhelpname="qhelp.dat",
             qname="q",
             auxname="a",
             plotq1name="plotq1"):
    """Generic code for plotting DoGPack output in matplotlib.

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
    from helper1 import read_params
    from helper1 import SampleBasis1
    from helper1 import read_qfile
    
    points_per_dir = int(points_per_dir_in)

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

    params = np.zeros(8, float) 
    GridType = read_params(outputdir,params)    
    meqn    = int(params[0])
    maux    = int(params[1])
    nplot   = int(params[2])
    meth1   = int(params[3])
    datafmt = int(params[4])
    mx      = int(params[5])
    xlow    = params[6]
    xhigh   = params[7]
    point_type = 1

    # TODO - this is hard-coded here
    meth1 = 1

    print("")
    print("        GridType = %s" % GridType        )
    print("  points_per_dir = %i" % points_per_dir  )
    print("      point_type = %i" % point_type      )
    print("       outputdir = %s" % outputdir       )
    print("       sorder    = %i" % meth1           )
    print("")

    # Grid information
    mx_old = mx
    mx = mx*points_per_dir
    dx_old = (xhigh-xlow)/float(mx_old)
    dx = (xhigh-xlow)/float(mx)

    xc = np.zeros(mx,float)
    xc[0] = xlow + 0.5*dx
    for i in range(1,mx):
        xc[i] = xc[i-1] + dx

    xc_old = np.zeros(mx_old,float)
    xc_old[0] = xlow + 0.5*dx_old
    for i in range(1,mx_old):
        xc_old[i] = xc_old[i-1] + dx_old

## ------------------------------------------------------------------------- ##
# Sample basis functions on mesh
# size of phi = (points_per_dir,meth1)
## ------------------------------------------------------------------------- ##
    phi = np.zeros((points_per_dir,meth1),float)
    SampleBasis1(points_per_dir,meth1,phi)

    q=-1;
    tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
    tmp2 = "".join((tmp1," ) ? "))
    m = raw_input(tmp2)
    print ""
    if (not m):
        m = 1
    else:
        m = int(m)

    if(m<1):
        print ""
        print "  Error, need m > 1,  m = ",m
        print ""
        exit(1)
    elif(m>meqn):
        print ""
        print "  Error, need m <=",meqn,",  m = ",m
        print ""
        exit(1)

    kn = 0;

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
        elif nf=="q":
            nf = -1
        else:
            nf = int(nf)
            n1 = nf

        if(n1>nplot):
            print ""
            print " End of plots "
            print ""
            n1 = nplot

        if(nf!=-1):

## ------------------------------------------------------------------------- ##
# Solution -- q
# solution should be found in file
#     outputdir/q[n1].dat
## ------------------------------------------------------------------------- ##
            qfile_tmp_tmp = "".join((str(n1+10000),".dat"))
            qfile_tmp = "q" + qfile_tmp_tmp[1:]
            qfile = "".join(("".join((outputdir,"/")),qfile_tmp))

            mtmp = mx_old*meqn*meth1
            qtmp = np.zeros(mtmp,float)   
            time = read_qfile(mtmp,qfile,qtmp)
            qcoeffs = np.reshape(qtmp,(meth1,meqn,mx_old))

            qsoln = np.zeros((mx,meqn),float);
            v1 = np.zeros(meth1,float)
            v2 = np.zeros(meth1,float)
            for i in range(1,mx_old+1):
                for me in range(1,meqn+1):
                    for ii in range(1,points_per_dir+1):
                        v1[:] = phi[ii-1,:]
                        v2[:] = qcoeffs[:,me-1,i-1]
                        tmp = 0.0
                        for k in range(0,meth1):
                            tmp = tmp + v1[k]*v2[k]
                        qsoln[(i-1)*points_per_dir+ii-1,me-1] = tmp

            if(maux>0):
                ## ------------------------------------------------------------------------- ##
                # Aux arrays -- aux
                # solution should be found in file
                #     outputdir/aux[n1].dat
                ## ------------------------------------------------------------------------- ##
                afile_tmp_tmp = "".join((str(n1+10000),".dat"))
                afile_tmp = auxname + afile_tmp_tmp[1:]
                afile = "".join(("".join((outputdir,"/")),afile_tmp))

                mtmp = mx_old*maux*meth1
                atmp = np.zeros(mtmp,np.float64)   
                time = read_qfile(mtmp,afile,atmp)
                acoeffs = np.reshape(atmp,(meth1,maux,mx_old))

                auxsoln = np.zeros((mx,maux),np.float64);
                v1 = np.zeros(meth1,np.float64)
                v2 = np.zeros(meth1,np.float64)
                for i in range(1,mx_old+1):
                    for me in range(1,maux+1):
                        for ii in range(1,points_per_dir+1):
                            v1[:] = phi[ii-1,:]
                            v2[:] = acoeffs[:,me-1,i-1]
                            tmp = 0.0
                            for k in range(0,meth1):
                                tmp = tmp + v1[k]*v2[k]
                            auxsoln[(i-1)*points_per_dir+ii-1,me-1] = tmp
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
    If executed at command line prompt, simply run plotdog1
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

#    This is how you could check if the user provided enough arguments...
#    if not opts.infile or not opts.outfile:
#        parser.error('Both options -i and -o are required. Try -h for help.')

### The old way of calling plotdog1: ###
#   import sys
#   args = sys.argv[1:]   # any command line arguments
#   plotdog1(*args)
