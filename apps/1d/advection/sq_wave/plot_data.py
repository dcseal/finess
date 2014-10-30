import numpy as np
import matplotlib.pyplot as plt

fname = 'data-good.dat'
data  = np.loadtxt( fname )
cfl_good   = data[:,0]
tvd_good   = data[:,1]

fname = 'data-bad.dat'
data  = np.loadtxt( fname )
cfl_bad   = data[:,0]
tvd_bad   = data[:,1]

fig = plt.figure()

ax1  = fig.add_subplot(111)
ax1.plot( cfl_good, tvd_good, 'bo' )
ax1.plot( cfl_bad, tvd_bad, 'kx' )
ax1.set_xlim( (min(cfl_good), max( cfl_good ) ) )
ax1.set_ylim( (-0.1, max( max( tvd_good ), max( tvd_bad ) )+0.1 ) )

ax1.set_xlabel('CFL number')
ax1.set_ylabel('Increase in Total Variation')

#ax2  = fig.add_subplot(221)
#ax2.plot( cfl_good, tvd_good, 'bo' )
#ax2.plot( cfl_bad, tvd_bad, 'kx' )
#ax2.set_xlim( (min(cfl_good), max( cfl_good ) ) )
#ax2.set_ylim( (-0.1, max( max( tvd_good ), max( tvd_bad ) )+0.1 ) )


fig.show()


