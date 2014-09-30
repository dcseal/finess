"""Provides functions for dealing with output of 2D apps.
"""

def read_qa(params, i_output):
    """Returns (t, q, aux) from (i_output)-th frame."""
    output_dir      = params["finess", "output_dir"]
    q_filename      = output_dir + '/' + ('q%.4d.dat' % i_output)
    aux_filename    = output_dir + '/' + ('a%.4d.dat' % i_output)
    
    mx      = params['grid', 'mx']
    meqn    = params['finess', 'meqn']
    maux    = params['finess', 'maux']

    from pylab import empty
    q = empty([mx, meqn])
    q_file = open(q_filename, 'r')
    t_q = float(q_file.readline())
    for m in range(meqn):
        for i in range(mx):
            q[i, m] = float(q_file.readline())
    q_file.close()
    
    aux         = empty([mx, maux])
    aux_file    = open(aux_filename, 'r')
    t_aux       = float(aux_file.readline())
    for m in range(maux):
        for i in range(mx):
            aux[i, m] = float(aux_file.readline())
    aux_file.close()
    
    assert t_q == t_aux, 'Inconsistent times from %s and %s' % (q_file, aux_file)
    
    return t_q, q, aux

