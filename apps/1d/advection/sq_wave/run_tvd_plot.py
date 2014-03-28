#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import numpy as np

def main( ):
    '''Open folders output00[nums] and extract the total variation from these
    files.
'''

    # parser for the parameters.ini files (needed to read in the cfl number)
    import ConfigParser, os
    config = ConfigParser.ConfigParser()

    my_dictionary = {}
    i = 0
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        folder = (os.getcwd() + '/output%(dir_num)04i/') % my_dictionary
        print folder
        if( not os.path.exists(folder) ):
#           print('    Did not find folder: %s' % folder )
            break

        my_dictionary['curr_folder'] = folder
        directory_num = i

        # Grab the current cfl number:
        config.readfp(open(folder + '/parameters.ini'))
        cfl_now = float( config.get('dogParams','cflv(2)') ) )

        try:
            qex  = np.loadtxt(folder + "/q0000.dat")[1:]
            qapp = np.loadtxt(folder + "/q0001.dat")[1:]
        except IOError:
            print('''Did not find the data file.
Please Wait for simulation to finish running.''')
            break

        i = i + 1

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
