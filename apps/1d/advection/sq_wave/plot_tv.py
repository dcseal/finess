import numpy as np
import matplotlib.pyplot as plt


fname = 'output/total-variation.dat'
fid   = open( fname, 'r' )

# read in the data
data = np.loadtxt( fname )
t  = data[:,0]
tv = data[:,1]

# plot the data
plt.figure()
plt.plot( t, tv, 'bo' )
plt.plot( t, 2.0*np.ones(t.shape), 'k--' )

# final call for all python plotting stuff
plt.show()
