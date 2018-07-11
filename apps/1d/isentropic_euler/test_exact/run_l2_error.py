#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np
import csv,itertools

def main( ):
    '''Write some help documentation here
'''

    print "# leading comments can be given a '#' character"
    my_dictionary = {}
    old_err = i = 0
    error=[['mx','L2','order']]
    o1=5
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        folder_RK = (os.getcwd() + '/output_Runge-Kutta_%(dir_num)03i/')  % my_dictionary 
        folder_DT = (os.getcwd() + '/output_User-Defined_%(dir_num)03i/') % my_dictionary 

        if( not os.path.exists(folder_RK) ):
            print 'did not find folder: %s' % folder
            break

        my_dictionary['curr_folder'] = folder_RK

        # we want to do:
        #   data = open('finess.data','w')
        #   print >> data, finess_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        directory_num = i

        try:
            qex  = np.loadtxt(folder_RK + "/q0001.dat")[1:]
            qapp = np.loadtxt(folder_DT + "/q0001.dat")[1:]
        except IOError:
            print('''Did not find the data file. Please Wait for simulation to finish running.''')
            break
 
        diff = qex - qapp

        new_err = np.sqrt( np.dot(diff,diff) ) / np.sqrt( np.dot(qex,qex) )

        r1 = '%(mx).f   L2-error = %(new).3e;  ' % {'old': old_err, 'new' : new_err, 'mx' : 0.5*qex.shape[0]}
        error1=[]
        
        if( old_err > 0 and new_err > 0 ):
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err/new_err), 1.2) } 
            error1=[0.5*qex.shape[0],new_err,log(old_err/new_err,1.2)]
        else:
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 
            error1=[0.5*qex.shape[0],new_err,0.000]

        error=error+[error1]
        print result
        old_err = new_err
        thefile = open( 'file%d.csv'%o1 , 'w' )
        for list1 in itertools.islice(error,1,None):
           writer = csv.writer(thefile)
           writer.writerow([list1[0],list1[1],list1[2]])
           #thefile.write("%s\n" % list1)

        old_err = new_err

        i = i + 1

if __name__ == '__main__':
    """TODO - using optparse is apparently deprecated in favor of argparse.
    This should be updated and/or deleted.  Right now, this script doesn't take
    in input arguments anyway. -DS 8/1/17"""

    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
