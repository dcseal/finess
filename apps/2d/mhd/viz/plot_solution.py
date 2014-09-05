def read_params(output_dir = 'output'):
    import os
    qhelp_filename = os.path.join(output_dir, 'qhelp.dat')
    qhelp_file = open(qhelp_filename, 'r')
    qhelp_lines = qhelp_file.readlines()
    qhelp_file.close()


    import string
    qhelp_params = map(lambda s: string.split(s)[0], qhelp_lines)
    
    ndims = int(qhelp_params[0])
    mesh_type = qhelp_params[1]
    meqn = int(qhelp_params[2])
    maux = int(qhelp_params[3])
    nout = int(qhelp_params[4])
    mx = int(qhelp_params[7])
    my = int(qhelp_params[8])
    xlow = float(qhelp_params[9])
    xhigh = float(qhelp_params[10])
    ylow = float(qhelp_params[11])
    yhigh = float(qhelp_params[12])
    dx = (xhigh-xlow) / float(mx)
    dy = (yhigh-ylow) / float(my)

    params = dict()
    params['ndims'] = ndims
    params['mesh_type'] = mesh_type
    params['meqn'] = meqn
    params['maux'] = maux
    params['nout'] = nout
    params['mx'] = mx
    params['my'] = my
    params['xlow'] = xlow
    params['xhigh'] = xhigh
    params['ylow'] = ylow
    params['yhigh'] = yhigh
    params['dx'] = dx
    params['dy'] = dy

    return params


def finess_2d_meshgrid(params):
    mx = params['mx']
    my = params['my']
    xlow = params['xlow']
    xhigh = params['xhigh']
    ylow = params['ylow']
    yhigh = params['yhigh']
    dx = params['dx']
    dy = params['dy']
    
    from pylab import meshgrid, linspace, empty
    X, Y = meshgrid(linspace(xlow + 0.5*dx, xhigh - 0.5*dx, mx), linspace(ylow + 0.5*dy, yhigh - 0.5*dy, my))
    return X, Y


def read_qa(params, i_output, output_dir = 'output'):
    """Returns (t, q, aux) from (i_output)-th frame."""
    from numpy import empty

    q_filename = output_dir + '/' + ('q%.4d.dat' % i_output)
    aux_filename = output_dir + '/' + ('a%.4d.dat' % i_output)
    
    mx = params['mx']
    my = params['my']
    meqn = params['meqn']

    q = empty([mx, my, meqn])
    q_file = open(q_filename, 'r')
    q_iter = iter(map(float, q_file))
    q_file.close()
    t_q = q_iter.next()
    for m in range(meqn):
        for j in range(my):
            for i in range(mx):
                q[i, j, m] = q_iter.next()
    

    maux = params['maux']
    aux = empty([mx, my, maux])
    aux_file = open(aux_filename, 'r')
    aux_iter = iter(map(float, aux_file))
    aux_file.close()
    t_aux = aux_iter.next()
    assert t_q == t_aux, 'Inconsistent times from %s and %s' % (q_file, aux_file)
    for m in range(maux):
        for j in range(my):
            for i in range(mx):
                aux[i, j, m] = aux_iter.next() 
    
    
    return t_q, q, aux

    
params = read_params()

X, Y = finess_2d_meshgrid(params)

t, q, aux = read_qa(params, 10)

from pylab import figure, show
from mpl_toolkits.mplot3d import axes3d, Axes3D            
fig = figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, q[:,:,0], rstride=2, cstride=2)
show()

