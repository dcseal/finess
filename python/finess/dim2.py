def meshgrid(params):
    """Returns meshgrid (a pair (X, Y)) that can be used for plotting."""
    assert params['finess', 'ndims'] == 2
    mx = params['grid', 'mx']
    my = params['grid', 'my']
    xlow = params['grid', 'xlow']
    xhigh = params['grid', 'xhigh']
    ylow = params['grid', 'ylow']
    yhigh = params['grid', 'yhigh']
    
    dx = (xhigh-xlow) / float(mx)
    dy = (yhigh-ylow) / float(my)
    
    from pylab import meshgrid, linspace
    X, Y = meshgrid(linspace(xlow + 0.5*dx, xhigh - 0.5*dx, mx), linspace(ylow + 0.5*dy, yhigh - 0.5*dy, my))
    return X, Y



def read_qa(params, i_output, output_dir = 'output'):
    """Returns (t, q, aux) from (i_output)-th frame."""
    q_filename = output_dir + '/' + ('q%.4d.dat' % i_output)
    aux_filename = output_dir + '/' + ('a%.4d.dat' % i_output)
    
    mx = params['grid', 'mx']
    my = params['grid', 'my']
    meqn = params['finess', 'meqn']
    maux = params['finess', 'maux']

    from pylab import empty
    q = empty([mx, my, meqn])
    q_file = open(q_filename, 'r')
    t_q = float(q_file.readline())
    for m in range(meqn):
        for j in range(my):
            for i in range(mx):
                q[i, j, m] = float(q_file.readline())
    q_file.close()
    
    aux = empty([mx, my, maux])
    aux_file = open(aux_filename, 'r')
    t_aux = float(aux_file.readline())
    for m in range(maux):
        for j in range(my):
            for i in range(mx):
                aux[i, j, m] = float(aux_file.readline())
    aux_file.close()
    
    assert t_q == t_aux, 'Inconsistent times from %s and %s' % (q_file, aux_file)
    
    return t_q, q, aux
:wa

