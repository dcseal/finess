import numpy as np
import matplotlib.pyplot as plt

fname = 'output/total-variation.dat'
fid   = open( fname, 'r' )

# read in the data
data = np.loadtxt( fname )
t  = data[:,0]
tv = data[:,1]

print("Largest increse in total variation = %2.5e" % ( max( tv[1:] - tv[0:len(tv)-1] ) ) )
pltsoln = True

print('tv[0] = ', tv[0])

if( pltsoln ):
    # plot the data
    plt.figure()
    plt.plot( t, tv-tv[0], 'bo' )
    plt.plot( t, np.zeros(t.shape), 'k--' )
#   plt.ylim( min( min(tv), 1.5 ), max( max(tv), 2.8 ) )
#   plt.ylim( min( min(tv), 1.5 ), 2.8 )

    # final call for all python plotting stuff
    plt.show()
