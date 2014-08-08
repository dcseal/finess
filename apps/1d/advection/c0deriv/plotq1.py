#----------------------------------------------------------
def plotq1(m, meth1, meqn, mx, time, xc, qsoln, auxsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # Exact solution:
    xold  = ( (1. + xc - time)%2. ) - 1.
    qex   = 0.*qsoln[:,0]   # Exact solution
    qexp  = 0.*qsoln[:,0]   # Derivative of exact solution
    for (i,x) in enumerate(xold):
        if( 0.25 < x and x < 0.4 ):
            qex [i] = (x-0.25)/0.075
            qexp[i] = 1./0.075
        elif( 0.4 <= x and x < 0.6 ):
            qex[i] = 2.0
        elif( 0.6 <= x and x < 0.75 ):
            qex[i]  = (0.75-x)/0.075
            qexp[i] = -1./0.075

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.1, 2.1])
    plt.plot(xc,qsoln[:,m],'bo')
    plt.plot(xc,qex, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    # Derivative of solution
    _stencil = [-2,-1,0,+1,+2]  # 5-point stencil
    _mbc     = 3                # required number of ghost-cells

    dx   = float( (xc[1] - xc[0]) )
    mbc  = 2
    qbig = np.concatenate( (qsoln[mx-mbc:,0], qsoln[:,0], qsoln[:mbc,0] ) )

    # Six point stencil (one 5 point for left, and one 5 point for right)
    Im2 = range(0, (mx+mbc-2))
    Im1 = range(1, (mx+mbc-1))
    I   = range(2, (mx+mbc  ))
    Ip1 = range(3, (mx+mbc+1))
    Ip2 = range(4, (mx+mbc+2))

    qim2 = qbig[Im2]
    qim1 = qbig[Im1]
    qi   = qbig[I  ]
    qip1 = qbig[Ip1]
    qip2 = qbig[Ip2]

    diffq = 1./(12.*dx)*( qim2-8.*( qim1 - qip1 ) - qip2 )

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim( [xc[0], xc[mx-1]] )
    plt.gca().set_ylim( [-15., 15.] )
    print( xc.shape, diffq.shape )
    plt.plot(xc, diffq, 'bo')
    plt.plot(xc,qexp, '-r', linewidth=2.0)
    tmp1 = "".join(("dq/dx(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()


#----------------------------------------------------------
