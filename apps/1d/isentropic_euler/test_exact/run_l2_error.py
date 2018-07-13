#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np
import csv,itertools

def main( ):
    '''Compute relative L2-errors for a numerical simulation and print the
    results to a csv file.
'''

    print("# leading comments can be given a '#' character")
    my_dictionary = {}
    error=[['mx','L2','order']]

    for dt_order in [5,7,9]:

        my_dictionary['dt_order'] = dt_order

        old_err = i = 0
        while( 1 ):

            directory_num = my_dictionary['dir_num'] =  i

            folder_RK = (os.getcwd() + '/output_Runge-Kutta_%(dir_num)03i/')  % my_dictionary 
            folder_DT = (os.getcwd() + '/output_User-Defined_Order%(dt_order)02i_%(dir_num)03i/') % my_dictionary 
            print(folder_DT)

            if( not os.path.exists(folder_RK) ):
                print('did not find folder: %s' % folder_RK)
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
                print('Did not find the data file. Please Wait for simulation to finish running.')
                break
     
            diff = qex - qapp

            new_err = np.sqrt( np.dot(diff,diff) ) / np.sqrt( np.dot(qex,qex) )

            mx_now = int( 0.5*qex.shape[0] )
            r1 = '%(mx).i   L2-error = %(new).3e;  ' % {'old': old_err, 'new' : new_err, 'mx' : mx_now }
            error1=[]
            if( old_err > 0 and new_err > 0 ):
                print('old_err, new_err, mx_now, mx_old = ', old_err, new_err, mx_now, mx_old )
                bs = float(mx_now)/float(mx_old)
                result = r1 + '   log2(ratio) = %(rat).3f' % {'rat' : log( (old_err/new_err), bs ) }
#               error1 = [ mx_now, new_err, log(old_err/new_err, bs ) ]
                error1 = [ mx_now, new_err]

                with closing(open('file%i.dat'%dt_order,'a')) as data:
                    print >> data, '%d %2.15e' % tuple(error1)

            else:
                result = r1 + '   log2(ratio) = %(rat).3f' % {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 
                error1 = [ mx_now, new_err, 0.000 ]
                error1 = [ mx_now, new_err]

                with closing(open('file%i.dat'%dt_order,'w')) as data:
                    print >> data, '%d %2.15e' % tuple(error1)

            print(result)
            old_err = new_err
            mx_old = mx_now

#           error=error+[error1]
#           thefile = open( 'file%d.csv' % dt_order, 'w' )
#           for list1 in itertools.islice(error,1,None):
#              writer = csv.writer(thefile)
#              writer.writerow([l for l in list1])


            old_err = new_err

            i = i + 1

if __name__ == '__main__':
    """TODO - implement optional parsing."""

    main( )
