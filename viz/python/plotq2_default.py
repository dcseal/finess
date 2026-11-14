#----------------------------------------------------------
def plotq2(xc, yc, xl, yl, qsoln, auxsoln, time, ini_params, m ):
    """Default 2D plotting file.

    This routine provides a handle that a user can swap out in a single
    application by placing a file called plotq2.py in the application
    directory.  Here, we provide a minimal working example that allows one
    to plot results.

    xc, yc    - cell centers, shape (mx, my)
    xl, yl    - cell corners, shape (mx+1, my+1)
    qsoln     - solution, shape (mx, my, meqn)
    auxsoln   - aux arrays, shape (mx, my, maux), or 0.0 when maux == 0
    m         - (zero-based) component of q to plot
    """

    import matplotlib.pyplot as plt

    mx = ini_params['mx']
    my = ini_params['my']

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim(xl[0,0],xl[mx,0])
    plt.gca().set_ylim(yl[0,0],yl[0,my])
    plt.pcolor(xl,yl,qsoln[:,:,m])
    plt.colorbar()
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [FINESS]"))
    plt.title(title)
    plt.draw()

#----------------------------------------------------------
