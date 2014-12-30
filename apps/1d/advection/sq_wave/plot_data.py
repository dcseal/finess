import numpy as np
import matplotlib.pyplot as plt

fname         = 'data1.dat'
data          = np.loadtxt( fname )
cfl1_vec      = data[:,0]
tvd1_vec      = data[:,1]

fname         = 'data2.dat'
data          = np.loadtxt( fname )
cfl2_vec      = data[:,0]
tvd2_vec      = data[:,1]

fig = plt.figure()

ax1  = fig.add_subplot(111)
#   ax1.plot( cfl_good, tvd_good, 'bo' )
#   ax1.plot( cfl_bad, tvd_bad, 'kx' )
#   ax1.set_xlim( (min(cfl_good), max( cfl_good ) ) )
#   ax1.set_ylim( (-0.1, max( max( tvd_good ), max( tvd_bad ) )+0.1 ) )

#   ax1.set_xlabel('CFL number')
#   ax1.set_ylabel('Increase in Total Variation')

#ax2  = fig.add_subplot(221)
#ax2.plot( cfl_good, tvd_good, 'bo' )
#ax2.plot( cfl_bad, tvd_bad, 'kx' )
#ax2.set_xlim( (min(cfl_good), max( cfl_good ) ) )
#ax2.set_ylim( (-0.1, max( max( tvd_good ), max( tvd_bad ) )+0.1 ) )

ax1.semilogy( cfl1_vec, tvd1_vec, 'bo' )
ax1.set_xlim( (min(cfl1_vec), max( cfl1_vec ) ) )
ax1.set_ylim( (-0.1*max( abs(tvd1_vec) ), 1.1*max( tvd1_vec )) ) 

ax1.set_xlabel('CFL number')
ax1.set_ylabel('Increase in Total Variation')
ax1.set_title('Square wave advection')

ax1.semilogy( cfl2_vec, tvd2_vec, 'g+' )
ax1.set_xlim( (min(cfl2_vec), max( cfl2_vec ) ) )
ax1.set_ylim( (-0.1*max( abs(tvd2_vec) ), 1.1*max( tvd2_vec )) ) 

ax1.legend( ('No Shu-Osher', 'Shu-Osher'), loc=0 )

fig.show()

fig.savefig('Comparison of two methods.pdf')
