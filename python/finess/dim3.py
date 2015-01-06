"""This module provides functions for dealing with output of 3D apps.
"""




def read_qa(params, i_output):
    """Returns (t, q, aux) from (i_output)-th frame."""
    datafmt = params["finess", "datafmt"]
    if datafmt == "ASCII":
        return _read_qa_ascii(params, i_output)
    elif datafmt == "Silo":
        return _read_qa_silo(params, i_output)
    else:
        raise RuntimeError("datafmt %s not supported" % datafmt)


def _read_qa_ascii(params, i_output):
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



def _read_qa_silo(params, i_output):
    """Silo reader."""
    output_dir = params["finess", "output_dir"]
    qa_filename = output_dir + '/' + ('qa%.4d.silo' % i_output)
    mx = params['grid', 'mx']
    my = params['grid', 'my']
    mz = params['grid', 'mz']
    meqn = params['finess', 'meqn']
    maux = params['finess', 'maux']

    from pylab import empty
    q = empty([mx, my, mz, meqn])
    aux = empty([mx, my, mz, maux])
   
    from pyvisfile import silo
    qafile = silo.SiloFile(qa_filename, create = False)

    quadmesh = qafile.get_quadmesh("quadmesh")
    t = quadmesh.dtime
    
    from pylab import empty
    q = empty([mx, my, mz, meqn])
    aux = empty([mx, my, mz, maux])
    
    for i in range(meqn):
        q[:, :, :, i] = qafile.get_quadvar("q%d" % (i+1)).vals[0]
    if maux > 0:
        for i in range(maux):
            aux[:, :, :, i] = qafile.get_quadvar("a%d" % (i+1)).vals[0]
        
    return t, q, aux

