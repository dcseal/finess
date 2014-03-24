import numpy as np
import matplotlib.pyplot as plt

fname = 'output/total-variation.dat'
fid   = open( fname, 'r' )

# read in the data
data = np.loadtxt( fname )
t  = data[:,0]
tv = data[:,1]

print("Largest increse in total variation = %2.5e" % ( max( tv[1:] - tv[0:len(tv)-1] ) ) )
pltsoln = False

if( pltsoln ):
    # plot the data
    plt.figure()
    plt.plot( t, tv, 'bo' )
    plt.plot( t, 2.0*np.ones(t.shape), 'k--' )
#   plt.ylim( min( min(tv), 1.5 ), max( max(tv), 2.8 ) )
    plt.ylim( min( min(tv), 1.5 ), 2.8 )

    # final call for all python plotting stuff
    plt.show()
