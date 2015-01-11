# Generate exact solutions for Burgers' equation


def read_q(parameters_ini_filename):
    from finess.params.util import read_params
    from finess.dim3 import read_qa
    from generate_iniparams import parameter_list
    params = read_params(parameters_ini_filename, parameter_list)
    nout = params['finess', 'nout']
    t, q, aux = read_qa(params, nout)
    return q

   

def log2_adjacent_ratio(error_list):
    order_list = []
    from numpy import log2
    for i in range(len(error_list) - 1):
        order_list.append(log2(error_list[i] / error_list[i+1]))
    return order_list


parameters_ini_filename_list = \
    ["parameters%02d.ini" % i for i in [0, 1, 2]]


error_list = []

from numpy import abs, max, sum
from numpy import sin, cos, pi
for parameters_ini_filename in parameters_ini_filename_list:
    from numpy import empty
    from finess.params.util import read_params
    from generate_iniparams import parameter_list
    params = read_params(parameters_ini_filename, parameter_list)
    print "Processing %s" % parameters_ini_filename
    q_original = read_q(parameters_ini_filename)
    q = q_original 
    from finess.viz.dim3 import meshgrid
    mx = params["grid", "mx"]
    my = params["grid", "my"]
    mz = params["grid", "mz"]
    t = params["finess", "tfinal"]
    X, Y, Z = meshgrid(params)
    qexact = empty((mx, my, mz, 8))
    qexact[:, :, :, 1 - 1] = 1.0
    qexact[:, :, :, 3 - 1] = 0.0
    qexact[:, :, :, 4 - 1] = 0.1 * sin(2.0*pi * (Y+t))
    qexact[:, :, :, 2 - 1] = 0.1 * cos(2.0*pi * (Y+t))
    qexact[:, :, :, 7 - 1] = 1.0
    qexact[:, :, :, 8 - 1] = 0.1 * sin(2.0*pi * (Y+t))
    qexact[:, :, :, 6 - 1] = 0.1 * cos(2.0*pi * (Y+t))
    pexact = 0.1
    # We don't care about energy, because we don't test the accuracy
    # of that.
    # qexact[:, :, :, 5 - 1]
    L1_error = sum(abs(q[:,:,:, 2 - 1]/q[:,:,:, 1 - 1] - qexact[:,:,:, 2 - 1]/qexact[:,:,:, 1 - 1])) / sum(abs(qexact[:,:,:, 2 - 1]/qexact[:,:,:, 1 - 1]))   #u1
    #L1_error = sum(abs(q[:,:,:, 4 - 1]/q[:,:,:, 1 - 1] - qexact[:,:,:, 4 - 1]/qexact[:,:,:, 1 - 1])) / sum(abs(qexact[:,:,:, 4 - 1]/qexact[:,:,:, 1 - 1]))   #u3
    #L1_error = sum(abs(q[:,:,:, 6 - 1] - qexact[:,:,:, 6 - 1])) / sum(abs(qexact[:,:,:, 6 - 1]))  #B1
    #L1_error = sum(abs(q[:,:,:, 8 - 1] - qexact[:,:,:, 8 - 1])) / sum(abs(qexact[:,:,:, 8 - 1]))  #B3
    error_list.append(L1_error)
    #error = abs(q[:,:,:,0]-qexact[:,:,:,0])
    #error_list.append(max(error))

order_list = log2_adjacent_ratio(error_list)

print "L1, error and order. u1"
print error_list
print order_list
