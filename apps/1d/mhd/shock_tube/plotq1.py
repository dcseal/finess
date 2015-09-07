#----------------------------------------------------------
def plotq1(m, meth1, meqn, mx, time, xc, qsoln, auxsoln):
    """Default plotting file.

    This routine provides a handle that a user can swap out in a single
    application.  Here, we provide a minimal working example that allows one
    to plot results.  For most problems, a user will likely want to swap this
    out for a tailored plotting script.
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    from math import fabs

    # Indices for components
    ind_rho       = 0; #  rho
    ind_M1        = 1; #  1-momentum
    ind_M2        = 2; #  2-momentum
    ind_M3        = 3; #  3-momentum
    ind_N         = 4; #  energy
    ind_B1        = 5; #  1-magnetic field
    ind_B2        = 6; #  2-magnetic field
    ind_B3        = 7; #  3-magnetic field
    ind_entropy_e = 8; #  entropy tracking

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
#   plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    plt.plot(xc,qsoln[:,ind_rho],'b-')
    tmp1 = "".join(("Density) at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [FINESS]"))
    plt.title(title)
    plt.draw()
    plt.savefig('MHD-Density.pdf')

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
#   plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    plt.plot(xc,qsoln[:,ind_B2],'b-')
    tmp1 = "".join(("Magnetic Field) at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [FINESS]"))
    plt.title(title)
    plt.draw()
    plt.savefig('MHD-Magnetic-Field.pdf')

#----------------------------------------------------------
