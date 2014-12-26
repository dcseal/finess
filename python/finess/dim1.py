"""This module provides functions for dealing with output of 1D apps.
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
    """Plain text reader."""
    output_dir = params["finess", "output_dir"]
    q_filename = output_dir + '/' + ('q%.4d.dat' % i_output)
    aux_filename = output_dir + '/' + ('a%.4d.dat' % i_output)
    
    mx = params['grid', 'mx']
    meqn = params['finess', 'meqn']
    maux = params['finess', 'maux']

    from pylab import empty
    q = empty([mx, meqn])
    q_file = open(q_filename, 'r')
    t_q = float(q_file.readline())
    for m in range(meqn):
        for i in range(mx):
            q[i, m] = float(q_file.readline())
    q_file.close()
    
    aux = empty([mx, maux])
    aux_file = open(aux_filename, 'r')
    t_aux = float(aux_file.readline())
    for m in range(maux):
        for i in range(mx):
            aux[i, m] = float(aux_file.readline())
    aux_file.close()
    
    assert t_q == t_aux, 'Inconsistent times from %s and %s' % (q_file, aux_file)
    
    return t_q, q, aux


def _read_qa_silo(params, i_output):
    """Silo reader."""
    output_dir = params["finess", "output_dir"]
    qa_filename = output_dir + '/' + ('qa%.4d.silo' % i_output)
    mx = params['grid', 'mx']
    meqn = params['finess', 'meqn']
    maux = params['finess', 'maux']

    from pylab import empty
    q = empty([mx, meqn])
    aux = empty([mx, maux])
   
    from pyvisfile import silo
    qafile = silo.SiloFile(qa_filename, create = False)

    quadmesh = qafile.get_quadmesh("quadmesh")
    t = quadmesh.dtime
    
    from pylab import empty
    q = empty([mx, meqn])
    aux = empty([mx, maux])
    
    q_raw = qafile.get_quadvar("q")
    for i in range(meqn):
        q[:, i] = q_raw.vals[i]
    if maux > 0:
        aux_raw = qafile.get_quadvar("a")
        for i in range(maux):
            aux[:, i] = aux_raw.vals[i]
        
    return t, q, aux
