
#----------------------------------------------------------
def plotq2_cart(outputdir,nframe,
                m,meth1,meqn,time,
                points_per_dir,LegVals,
                xlow,xhigh,ylow,yhigh,
                mx,my,
                dx,dy,
                mx_old,my_old,
                dx_old,dy_old,
                xc,yc,
                xl,yl,
                qaug):
    
    import matplotlib.pyplot as plt
    import numpy as np
    from math import sqrt
    from math import pow
    from math import cos
    from math import sin
    from math import pi

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim(xl[0,0],xl[mx,0])
    plt.gca().set_ylim(yl[0,0],yl[0,my])
    plt.pcolor(xl,yl,qaug[:,:,m],vmin=0.0,vmax=1.0)
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    x0 = -0.1*cos(2.0*pi*time) + 0.50
    y0 =  0.1*sin(2.0*pi*time) + 0.50
    r = np.zeros(mx*my,float)
    qscat = np.zeros((mx*my,meqn),float)
    k = 0
    for i in range(0,mx):
        for j in range(0,my):            
            r[k] = sqrt(pow(xc[i,j]-x0,2)+pow(yc[i,j]-y0,2))
            for m in range(0,meqn):
                qscat[k,m] = qaug[i,j,m]
            k = k+1
    qscat_ex = np.zeros((mx*my,meqn),float)
    qex(mx,my,meqn,r,qscat_ex)
    
    ind = r.argsort()
    
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([0.0,0.5])
    plt.gca().set_ylim([0.0,1.0])
    plt.plot(r[ind],qscat_ex[ind,m],'k-')
    plt.plot(r[ind],qscat[ind,m],'bo')
    tmp1 = "".join(("Scattor plot of q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------


#----------------------------------------------------------
def qex(mx,my,meqn,r,qscat_ex):

    from math import pow
    from math import cos
    from math import pi
    
    for i in range(mx*my):
        if (r[i]<0.3):
            qscat_ex[i,0] = pow( cos(5.0/3.0*pi*r[i]) ,6)

#----------------------------------------------------------
