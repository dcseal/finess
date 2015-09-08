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

    qlow  = min(qsoln[:,m])
    qhigh = max(qsoln[:,m])
    if fabs(qhigh-qlow)>1.0e-7:            
        qeps  = 0.1*(qhigh-qlow)
    else:
        qeps = 0.1            

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    plt.plot(xc,qsoln[:,m],'b-')
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [FINESS]"))
    plt.title(title)
    plt.draw()
    plt.pause(0.1)
    
#----------------------------------------------------------
