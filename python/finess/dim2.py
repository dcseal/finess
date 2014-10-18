"""This module provides functions for dealing with output of 2D apps.
"""

def read_qa(params, i_output):
    """Returns (t, q, aux) from (i_output)-th frame."""
    output_dir = params["finess", "output_dir"]
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

