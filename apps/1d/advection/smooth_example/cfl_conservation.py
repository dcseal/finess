#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import numpy as np

def main( ):
    '''Open folders output00[nums] and extract the total variation and conserved
    quatities from these files.

    See also: ConSoln.cpp.
'''

    # parser for the parameters.ini files (needed to read in the cfl number)
    import ConfigParser, os
    config = ConfigParser.ConfigParser()

    # Count the number of output directories:
    my_dictionary = {}
    i = 1
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        folder = (os.getcwd() + '/output%(dir_num)04i/') % my_dictionary
        if( not os.path.exists(folder) ):
            break
        i = i+1

    tvd_vec = np.zeros( (i,1) )
    print("found %i folders" % len(tvd_vec) )

    print("-----------------------------------------------------------------------")
    print(" CFL;        Mass;          L1-norm;        L2-norm;       Linf-norm   ")
    print("-----------------------------------------------------------------------")

    for i in range( len(tvd_vec) ):

        directory_num = my_dictionary['dir_num'] = i
        folder = (os.getcwd() + '/output%(dir_num)04i/') % my_dictionary

        # Grab the current cfl number:
        config.readfp(open(folder + '/parameters.ini'))
        cfl_now = float( config.get('dogParams','cflv(2)') )

        data   = np.loadtxt( folder + '/conservation.dat' )
        t      = data[:,0]
        mass   = data[:,1]
        l1norm = data[:,2]
        l2norm = data[:,3]
        linorm = data[:,4]

        dm   = max( mass   ) - min( mass   )
        d1   = max( l1norm ) - min( l1norm )
        d2   = max( l2norm ) - min( l2norm )
        di   = max( linorm ) - min( linorm )
        print("%2.3f;    %2.5e;    %2.5e;    %2.5e;    %2.5e" % ( cfl_now, dm, d1, d2, di ) )

if __name__ == '__main__':
    main.__doc__

    main( )
