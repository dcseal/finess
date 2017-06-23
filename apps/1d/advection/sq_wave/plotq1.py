#----------------------------------------------------------
def plotq1(xc, qsoln, auxsoln, time, ini_params, m ):
    
    import matplotlib.pyplot as plt
    import numpy as np

    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    nplot    = ini_params['nout']
    mx       = ini_params['mx']

    # Exact solution:
    width = 0.5
    xold = ( (xc+1.0 - time)%2. ) - 1.0
    qex  = 0.*qsoln[:,0]
    for (i,x) in enumerate(xold):
        if( abs(x) < width ):
            qex[i] = 1.0
#           qex[i] = np.cos(np.pi*(x-0.5)/width)**6

    print("Linf-error = %2.15e\n" % max( abs(qsoln[:,0]-qex[:]) ) )
    print('total mass = %f\n'     % ( (xc[1]-xc[0])*sum( qsoln[:,0] ) ) )
#   print('sum = %f\n' % ( ( xc[1]-xc[0] ) * sum(qsoln[:,0]) ) )

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
