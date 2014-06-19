#----------------------------------------------------------
def plotq1(m,meth1,meqn,mx,time,xc,qsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.1,1.1])
    plt.plot(xc,qsoln[:,m],'bo')
    plt.grid()
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    ##=======================================================================##
    # Add in the exact solution
    ##=======================================================================##
    qs_left_problem  = 0.1339745962155613
    qs_right_problem = 0.5

    # Derivative of flux function:
    def buckley_fp( q ):
        M = 1./3.
        if( time > 0. ):
            return (2.*M*q*(1.-q) ) / (q**2 + M*(1.-q)**2)**2
        else:
            return 0*q

    qsample = np.linspace(qs_right_problem, 1.0 )
    dt      = time

    ql = np.linspace(0., qs_left_problem )
    left_shock = -0.5 + dt * buckley_fp( qs_left_problem )
    fl = -0.5 + dt*buckley_fp( ql )

    # segment 1
    plt.plot( [-1, -0.5], [0, 0], '-r' )

    # segment 2
    plt.plot( fl, ql, '-r' );


    # segment 3
    plt.plot( fl, ql, '-r' );
    plt.plot( [left_shock, left_shock], [qs_left_problem, 1.0], '-r' );

    # segment 4
    plt.plot( [left_shock, 0], [1, 1], '-r' );


    # segment 5
    plt.plot( dt*buckley_fp(qsample), qsample, '-r'  );

    # segment 6
    plt.plot( [ dt*buckley_fp(qs_right_problem), dt*buckley_fp(qs_right_problem) ], 
          [ qs_right_problem, 0 ], '-r' );
    
    # segment 7
    plt.plot( [dt*buckley_fp(qs_right_problem), 1.0], [0, 0], '-r' );

    # segment 8
    plt.plot( [-1, 1], [0, 0], '--k' );

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([-0.08,0.02])
    plt.gca().set_ylim([0.9,1.1])
    plt.plot(xc,qsoln[:,m],'bo')
    plt.grid()
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    ##=======================================================================##
    # Add in the exact solution
    ##=======================================================================##
    qs_left_problem  = 0.1339745962155613
    qs_right_problem = 0.5

    # Derivative of flux function:
    def buckley_fp( q ):
        M = 1./3.
        if( time > 0. ):
            return (2.*M*q*(1.-q) ) / (q**2 + M*(1.-q)**2)**2
        else:
            return 0*q

    qsample = np.linspace(qs_right_problem, 1.0 )
    dt      = time

    ql = np.linspace(0., qs_left_problem )
    left_shock = -0.5 + dt * buckley_fp( qs_left_problem )
    fl = -0.5 + dt*buckley_fp( ql )

    # segment 1
    plt.plot( [-1, -0.5], [0, 0], '-r' )

    # segment 2
    plt.plot( fl, ql, '-r' );


    # segment 3
    plt.plot( fl, ql, '-r' );
    plt.plot( [left_shock, left_shock], [qs_left_problem, 1.0], '-r' );

    # segment 4
    plt.plot( [left_shock, 0], [1, 1], '-r' );


    # segment 5
    plt.plot( dt*buckley_fp(qsample), qsample, '-r'  );

    # segment 6
    plt.plot( [ dt*buckley_fp(qs_right_problem), dt*buckley_fp(qs_right_problem) ], 
          [ qs_right_problem, 0 ], '-r' );
    
    # segment 7
    plt.plot( [dt*buckley_fp(qs_right_problem), 1.0], [0, 0], '-r' );

    # segment 8
    plt.plot( [-1, 1], [0, 0], '--k' );


    data = [xc, qsoln[:,0]]
    np.savetxt( 'dg_picture.dat', np.transpose( data ) )
#----------------------------------------------------------
