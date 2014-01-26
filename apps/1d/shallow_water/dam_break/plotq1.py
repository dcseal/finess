from __future__ import print_function

#----------------------------------------------------------
def plotq1(m,meth1,meqn,mx,time,xc,qsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # HEIGHT
    m=0
    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([0.9,3.1])
    plt.plot(xc,qsoln[:,m],'bo')
    tmp1 = "".join(("Height at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    # VELOCITY
    m=1
    u = np.zeros(mx,float)
    for i in range(0,mx):
        u[i] = qsoln[i,1]/qsoln[i,0]
        
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.1,0.8])
    plt.plot(xc,u,'bo')
    tmp1 = "".join(("Velocity at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
    data = np.array( [xc, qsoln[:,0], qsoln[:,1] ] )
# For some reason, these calls have a bug in them:
#   print(data)
#   with open( 'dg_picture.dat', 'wb' ) as f:
#     print( '%.15e' % time, file=f )         # time instant on first row
#     np.savetxt( 'dg_picture.dat', np.transpose( data ) )
    np.savetxt( 'dg_picture.dat', np.transpose( data ) )
 

#----------------------------------------------------------
