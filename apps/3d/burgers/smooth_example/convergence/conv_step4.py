# Generate exact solutions for Burgers' equation


def read_q1(parameters_ini_filename):
    from finess.params.util import read_params
    from finess.dim3 import read_qa
    from generate_iniparams import parameter_list
    params = read_params(parameters_ini_filename, parameter_list)
    nout = params['finess', 'nout']
    t, q, aux = read_qa(params, nout)
    return q[:, :, :, 1 - 1]

def read_qexact(parameters_ini_file_name):
    from finess.params.util import read_params
    from generate_iniparams import parameter_list
    params = read_params(parameters_ini_filename, parameter_list)
    output_dir = params['finess', 'output_dir']
    import os
    qexact_filename = os.path.join(output_dir, "qexact.silo")
    from pyvisfile import silo
    with silo.SiloFile(qexact_filename, create=False) as qexact_file:
        qexact_raw = qexact_file.get_quadvar("q")
        qexact = qexact_raw.vals[0]
    return qexact

def log2_adjacent_ratio(error_list):
    order_list = []
    from numpy import log2
    for i in range(len(error_list) - 1):
        order_list.append(log2(error_list[i] / error_list[i+1]))
    return order_list


parameters_ini_filename_list = \
    ["parameters%02d.ini" % i for i in [0, 1, 2, 3]]


error_list = []

from numpy import abs, max, sum
for parameters_ini_filename in parameters_ini_filename_list:
    print "Processing %s" % parameters_ini_filename
    q = read_q1(parameters_ini_filename)
    qexact = read_qexact(parameters_ini_filename)
    L1_error = sum(abs(q - qexact)) / sum(abs(qexact))
    error_list.append(L1_error)
#    error = abs(q-qexact)
#    error_list.append(max(error))
order_list = log2_adjacent_ratio(error_list)

print "L1, error and order"
print error_list
print order_list
