"""This module provides functions for dealing with output of 3D apps.
"""

def read_qa(params, i_output):
    """Returns (t, q, aux) from (i_output)-th frame."""

    output_dir = params["finess", "output_dir"]
    q_filename = output_dir + '/' + ('q%.4d.dat' % i_output)
    aux_filename = output_dir + '/' + ('a%.4d.dat' % i_output)
   
    # Grid information
    mx = params['grid', 'mx']
    my = params['grid', 'my']
    mz = params['grid', 'mz']

    # Problem information
    meqn = params['finess', 'meqn']
    maux = params['finess', 'maux']

    from pylab import empty
    q = empty([mx, my, mz, meqn])
    q_file = open(q_filename, 'r')
    t_q = float(q_file.readline())
    for m in range(maux):
      for k in range(mz):
        for j in range(my):
          for i in range(mx):
                q[i, j, k, m] = float(q_file.readline())
    q_file.close()
    
    aux = empty([mx, my, mz, maux])
    aux_file = open(aux_filename, 'r')
    t_aux = float(aux_file.readline())
    for m in range(maux):
      for k in range(mz):
        for j in range(my):
          for i in range(mx):
            aux[i, j, k, m] = float(aux_file.readline())
    aux_file.close()
    
    assert t_q == t_aux, 'Inconsistent times from %s and %s' % (q_file, aux_file)
    
    return t_q, q, aux

