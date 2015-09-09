#-----------------------------------------------------------------------------#
def plotq1(xc, qsoln, auxsoln, time, ini_params, m ):
    """Default plotting file for Shallow water equations. 

    This routine provides a handle that a user can swap out in a single
    application.  Here, we provide a minimal working example that allows one
    to plot results.  For most problems, a user will likely want to swap this
    out for a tailored plotting script.
    """
    
    import matplotlib.pyplot as plt
    import numpy as np

    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    mx       = ini_params['mx']

    # HEIGHT
    m=0
    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([0.1,1.2])
    plt.plot(xc,qsoln[:,m],'bo')
    tmp1 = "".join(("Height at t = ",str(time)))
    title = "".join((tmp1,"     [FINESS]"))
    plt.title(title)
    plt.draw()

    # VELOCITY
    m=1
    u = np.zeros(mx,float)
    u[:] = qsoln[:,1] / qsoln[:,0]

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
#   plt.gca().set_ylim([-0.1,0.8])
    plt.plot(xc,u,'bo')
    tmp1 = "".join(("Velocity at t = ",str(time)))
    title = "".join((tmp1,"     [FINESS]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
