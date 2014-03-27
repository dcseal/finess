
#----------------------------------------------------------
def plotq1(m,meth1,meqn,mx,time,xc,qsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # Convenience accessors
    pi  = np.pi
    cos = np.cos
    sin = np.sin
    exp = np.exp

    # Exact solution:
    f      = np.exp( np.cos( 2.0*np.pi*xc ) )
    fx     = -2.*pi*exp(np.cos(2*pi*xc))*np.sin(2*pi*xc)
    fxx    = 4*pi**2*(-cos(2*pi*xc)**2 - cos(2*pi*xc) + 1)*exp(cos(2*pi*xc))

    # Derivatives of the initial conditions
    qex    = -fx
    qex_x  = -fxx

    qex_xx = 0. + xc
    for n in range( len(xc) ):
        if( xc[n] > 0.3 and xc[n] < 0.7 ):
            qex_xx[n] = 1.
        else:
            qex_xx[n] = 0.

    print("Linf-error = %2.15e" % max( abs(qsoln[:,0]-qex[:]) ) )

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
#   plt.gca().set_ylim([-0.1,1.1])
    plt.plot(xc,qsoln[:,0],'bo')
    plt.plot(xc,qex, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
#   plt.gca().set_ylim([-0.1,1.1])
    plt.plot(xc,qsoln[:,1],'bo')
    plt.plot(xc,qex_x, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)

    plt.figure(3)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
#   plt.gca().set_ylim([-0.1,1.1])
    plt.plot(xc,qsoln[:,2],'bo')
    plt.plot(xc,qex_xx, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
