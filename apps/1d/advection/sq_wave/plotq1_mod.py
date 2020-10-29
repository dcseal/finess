
#----------------------------------------------------------
def plotq1(xc, qsoln, auxsoln, time, ini_params, m ):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # (Pulled from plotfin1.py)
    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    nplot    = ini_params['nout']
    mx       = ini_params['mx'  ]

    # Exact solution:
    width = 0.5
    xold  = ( (1. + xc - time)%2. ) - 1.
    qex   = 0.*qsoln[:,0]
    for (i,x) in enumerate(xold):
        if( abs(x) < width ):
            qex[i] = 1.0

    print('total mass = %f\n' % (xc[1]-xc[0])*sum( qsoln[:,0] ) )

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.1,1.1])
    plt.plot(xc,qsoln[:,m],'bo')
    plt.plot(xc,qex, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
