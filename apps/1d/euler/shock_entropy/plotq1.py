from __future__ import print_function  # used for printing to file

#-----------------------------------------------------------------------------#
def plotq1(xc, qsoln, auxsoln, time, ini_params, m ):
    
    import matplotlib.pyplot as plt
    import numpy as np

    meqn     = ini_params['meqn']
    maux     = ini_params['maux']
    mx       = ini_params['mx']

    # Pull gas constant (or read it if it hasn't already been read)
    if( not ini_params.has_key('gamma') ):

        import ConfigParser
        config = ConfigParser.RawConfigParser()
        config.read( ini_params['output_dir'] + '/parameters.ini' )
        ini_params['gamma'] = config.getfloat('euler', 'gamma' )

    gamma   = ini_params['gamma']
    gp1     = gamma+1.
    gm1     = gamma-1.

    # DENSITY
    m=0
    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([0.39, 5.1])
    plt.plot(xc,qsoln[:,m],'bo')
    tmp1 = "".join(("Density at t = ",str(time)))
    title = "".join((tmp1,"     [FINESS]"))
    plt.title(title)
    plt.draw()

    # PRESSURE
    m=1

    press = (gamma-1)*(qsoln[:,4]-0.5*(
        qsoln[:,1]**2 + qsoln[:,2]**2 + qsoln[:,3]**2)/qsoln[:,0])
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.01, 12.10])
    plt.plot(xc,press,'bo')
    tmp1 = "".join(("Pressure at t = ",str(time)))
    title = "".join((tmp1,"     [FINESS]"))
    plt.title(title)
    plt.draw()
 
    # Velocity
    m=2
    u = np.zeros(mx,float)
    u = qsoln[:,1] / qsoln[:,0]
        
    plt.figure(3)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.11,  3.10])
    plt.plot(xc,u,'bo')
    tmp1 = "".join(("Velocity at t = ",str(time)))
    title = "".join((tmp1,"     [FINESS]"))
    plt.title(title)
    plt.draw()


    # print data to file.
    data = np.array( [xc, qsoln[:,0], qsoln[:,1], qsoln[:,4] ] )
    fmt   = '%1.15e'
    with open( 'soln.dat', 'wb' ) as f:
      #print( fmt % time, file=f )         # time instant on first row
      np.savetxt( 'soln.dat', np.transpose( data ), fmt=fmt )
 
#-----------------------------------------------------------------------------#
