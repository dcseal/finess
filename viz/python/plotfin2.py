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

    # Pull the personalized plotting routine, plotq1.py
    curr_dir = os.path.abspath("./")
    sys.path.append(curr_dir)
    plotq2_file  = os.path.abspath("plotq2.py")
    local_plotq1 = os.path.exists(plotq1_file)
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

            # USER SUPPLIED FUNCTION
#           plotq2_cart(outputdir,n1,
#                       m-1,meth1,meqn,time,
#                       points_per_dir,LegVals,
#                       xlow,xhigh,ylow,yhigh,
#                       mx,my,dx,dy,
#                       mx_old,my_old,dx_old,dy_old,
#                       xc,yc,xl,yl,qaug);
        
    plt.ioff()
    print ""

elif (GridType=="Unstructured"):
    plotq2_file  = os.path.abspath("plotq2_unst.py")
    local_plotq2 = os.path.exists(plotq2_file)
    if local_plotq2==False:
        from plotq2_unst_default import plotq2_unst
    else:
        from plotq2_unst import plotq2_unst

    print " Creating mesh ... "

    # READ-IN MESH INFO
    meshdir = "".join((outputdir,"/mesh_output"))
    TF = os.path.exists(meshdir)
    if TF==False:
        print ""
        print "    Directory not found, meshdir =",meshdir
        print ""
        return -1

    tnode = np.zeros((NumPhysElems,3),int)
    x  = np.zeros(NumPhysNodes,float)
    y  = np.zeros(NumPhysNodes,float)
    read_node(meshdir,NumPhysNodes,x,y)
    read_tnode(meshdir,NumPhysElems,tnode)

    xlow  = min(x)
    ylow  = min(y)
    xhigh = max(x)
    yhigh = max(y)

    tcounter = np.zeros(NumPhysNodes,int)
    set_tcounter(NumPhysElems,tnode,tcounter)

    # Add extra points and elements if points_per_dir>1
    p2 = points_per_dir*points_per_dir
    points_per_elem = ((points_per_dir+1)*(points_per_dir+2))/2
    zx = np.zeros(p2,float)
    zy = np.zeros(p2,float)
    
    if (points_per_dir>1):
        x_tmp = np.zeros(NumPhysElems*points_per_elem,float)
        y_tmp = np.zeros(NumPhysElems*points_per_elem,float)
        tnode_new = np.zeros((NumPhysElems*p2,3),int)
        NumPhysElems_new = 0
        NumPhysNodes_new = 0
        newsizes = np.zeros(2,int)            
        DivideUnstMesh(points_per_dir,NumPhysElems,x,y,tnode,
                       x_tmp,y_tmp,tnode_new,zx,zy,newsizes)
        NumPhysElems_new = newsizes[0]
        NumPhysNodes_new = newsizes[1]
        x_new = np.zeros(NumPhysNodes_new,float)
        y_new = np.zeros(NumPhysNodes_new,float)
        for ijk in range(0,NumPhysNodes_new):
            x_new[ijk] = x_tmp[ijk]
            y_new[ijk] = y_tmp[ijk]
        tcounter_new = np.zeros(NumPhysNodes_new,int)
        set_tcounter(NumPhysElems_new,tnode_new,tcounter_new)
    else:
        NumPhysElems_new = NumPhysElems
        NumPhysNodes_new = NumPhysNodes
        x_new = x
        y_new = y
        tnode_new = tnode
        tcounter_new = tcounter

    # Get physical midpoints of each element
    xmid = np.zeros(NumPhysElems_new,float)
    ymid = np.zeros(NumPhysElems_new,float)
    onethird = 1.0/3.0
    for i in range(0,NumPhysElems_new):
        xmid[i] = onethird*(x_new[tnode_new[i,0]]+x_new[tnode_new[i,1]]+x_new[tnode_new[i,2]])
        ymid[i] = onethird*(y_new[tnode_new[i,0]]+y_new[tnode_new[i,1]]+y_new[tnode_new[i,2]])

    # Sample Legendre polynomial on the midpoint of each element
    Mon2Leg = np.zeros((kmax,kmax),float)
    GetMonomialToLegendre(kmax,Mon2Leg)
    MonVals = np.zeros((kmax,p2),float)
    LegVals = np.zeros((kmax,p2),float)
    GetUnstLegendre(meth1,kmax,points_per_dir,zx,zy,Mon2Leg,MonVals,LegVals)

    print " Finished creating mesh. "
    print ""

    q=-1;
    tmp1 = "".join((" Which component of q do you want to plot ( 1 - ",str(meqn)))
    tmp2 = "".join((tmp1," ) ? "))
    m = raw_input(tmp2)
    print ""
    if (not m):
        m = 1
    else:
        m = int(m)

    if m<1:
        print ""
        print "  Error, need m > 1,  m = ",m
        print ""
        return -1
    elif m>meqn:
        print ""
        print "  Error, need m <=",meqn,",  m = ",m
        print ""
        return -1

    kn = 0;

    n = 0;
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

        if n1>nplot:
            print ""
            print " End of plots "
            print ""
            n1 = nplot

        if (nf!=-1):
            # Solution -- q
            # solution should be found in file
            #     outputdir/q[n1].dat

            qfile_tmp_tmp = "".join((str(n1+10000),".dat"))
            qfile_tmp = "q" + qfile_tmp_tmp[1:]
            qfile = "".join(("".join((outputdir,"/")),qfile_tmp))

            mtmp = NumElems*meqn*kmax
            qtmp = np.zeros(mtmp,float)   
            time = read_qfile(mtmp,qfile,qtmp)                
            qcoeffs_tmp = np.reshape(qtmp,(NumElems,meqn,kmax),'fortran')
            qcoeffs = qcoeffs_tmp[0:NumPhysElems,:,:]

            # Solution -- q
            qsoln_elem = np.zeros((NumPhysElems_new,meqn),float)
            qsoln      = np.zeros((NumPhysNodes_new,meqn),float)
            sample_state2_unst(NumPhysElems,p2,meqn,kmax,LegVals,qcoeffs,qsoln_elem)
            set_soln_at_node_values(NumPhysElems_new,NumPhysNodes_new,
                                    meqn,tcounter_new,tnode_new,qsoln_elem,qsoln)

            # USER SUPPLIED FUNCTION: Plotting function
            plotq2_unst(outputdir, n1, 
                        m-1, meqn, NumPhysElems_new, NumPhysNodes_new, 
                        xlow, xhigh, ylow, yhigh, time, 
                        x_new, y_new, tnode_new, 
                        qsoln, xmid, ymid, qsoln_elem)                
        
    plt.ioff()
    print ""
    
else:
    
    print ""
    print " Error in plotfin2.py: GridType = ",GridType," is not supported."
    print ""
    return -1
            


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
    args = parse_input( plotfin2.__doc__ )

    # Call the main 1D plotting routine
    plotfin2(args.outputdir, args.q_name, args.aux_name, args.plotq1_name )

