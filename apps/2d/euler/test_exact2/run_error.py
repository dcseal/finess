#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log,sqrt
import numpy as np

def main( ):
    '''Write some help documentation here
'''

    print "# leading comments can be given a '#' character"
    my_dictionary = {}
    old_err = i = 0
    old_err2 = 0
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        folder = (os.getcwd() + '/output_%(dir_num)03i/') % my_dictionary

#       print folder
        if( not os.path.exists(folder) ):
            print 'did not find folder: %s' % folder
            break

        my_dictionary['curr_folder'] = folder
        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        directory_num = i

        try:
            qex  = np.loadtxt(folder + "/q0000.dat")[1:]
            qapp = np.loadtxt(folder + "/q0001.dat")[1:]
        except IOError:
            print('''Did not find the data file.
Please Wait for simulation to finish running.''')
            break
 
        qlength = len(qex)/5
        m = sqrt(qlength)
        dx = dy = 10.0/m
        print 'm = %(mm)d' % {'mm':m}
        qex =  qex[:qlength]    # only density for this error
        qapp = qapp[:qlength]   # only density
        diff = qex - qapp

        new_err =  sum(abs(diff)) * dx * dy /100.0
        new_err2 = max(abs(diff)) # / max(abs(qex))

        r1 = 'L1-error = %(new).3e;  ' % {'old': old_err, 'new' : new_err}
        
        if( old_err > 0 and new_err > 0 ):
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err/new_err), 2) } 
        else:
            result = r1 + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 

        r2 = 'Linf-error = %(new).3e;  ' % {'old': old_err2, 'new' : new_err2}

        if( old_err2 > 0 and new_err2 > 0 ):
            result2 = r2 + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err2/new_err2), 2) } 
        else:
            result2 = r2 + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err2, 'new' : new_err2, 'rat' : (old_err2/new_err2) } 


 # This is exactly the format I want:
#{\normalsize $25$}   &   {\normalsize $1.747\times 10^{-4}$}   & {\normalsize --} &  {\normalsize $8.292\times 10^{-5}$} & {\normalsize --}  \\
       
        print result 
        print result2

        old_err = new_err
        old_err2 = new_err2

        i = i + 1

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
