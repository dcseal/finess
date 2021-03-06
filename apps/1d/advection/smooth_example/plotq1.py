#----------------------------------------------------------
def plotq1(xc, qsoln, auxsoln, time, ini_params, m ):
    
    import matplotlib.pyplot as plt
    import numpy as np

    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    mx       = ini_params['mx']

    # Exact solution:
    width = 2.0*0.25
    xold = (xc - time)%1.
    qex  = 0.*qsoln[:,0]
    for (i,x) in enumerate(xold):
        if( abs(x-0.5) < 0.5*width ):
            qex[i] = np.cos(np.pi*(x-0.5)/width)**6

    print("Linf-error = %2.15e" % max( abs(qsoln[:,0]-qex[:]) ) )

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.1,1.1])
    plt.plot(xc,qsoln[:,m],'bo')
    plt.plot(xc,qex, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [FINESS]"))
    plt.title(title)
    plt.draw()

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.003,0.003])
    plt.plot(xc,qsoln[:,m]-qex,'bo')
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [FINESS]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
