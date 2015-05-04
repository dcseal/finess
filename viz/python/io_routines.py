# The purpose of this module is to contain common IO (Input/Output) routines
# used for plotting unkowns in FINESS.
#
# Currently only supports 1D plotting.  
# See $(FINESS)/lib/[1-3]d/Output.cpp to see the order in which variables are saved to file.

#----------------------------------------------------------
def read_params(outputdir):
    """Routine used for parsing any qhelp.dat file."""

    import string

    params  = {}
    Fparams = "".join((outputdir,"/qhelp.dat"     ))
    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist   = string.split(linestring)
    ndims      = int(linelist[0])

    if(not ndims==1 ):
        print("\n Incorrect dimension, ndims must be 1. ndims = \n", ndims )
        return -1

    for k in range(6):
        linestring = Rparams.readline()
        linelist   = string.split(linestring)
        params[linelist[2]] = float(linelist[0])

    Rparams.close()
    return params

#----------------------------------------------------------
def read_qfile(mx, meqn, qfile ):
    """Read solution from file.  The format of this file is identical to that
    used by FINESS for output.
    """

    import string
    import numpy as np

    # solution to be returned
    qsoln = np.zeros((mx, meqn), float);

    # open file
    Rqfile = open(qfile,'r')

    # get time
    linestring = Rqfile.readline()
    linelist   = string.split(linestring)
    time       = float(linelist[0])
    
    # extract all point values
    k = 0
    for i in range(mx):
        for me in range(meqn):

            linestring  = Rqfile.readline()
            linelist    = string.split(linestring)
            qsoln[i,me] = float(linelist[0])

            k = k+1

    # close file
    Rqfile.close()

    # return time
    return time, qsoln
#----------------------------------------------------------
