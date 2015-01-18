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
    ["output%(i)02d/parameters%(i)02d.ini.dump" % {"i": i} for i in [0, 1, 2, 3, 4]]


error_list = []

from numpy import abs, max, sum, zeros, ones
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
    theta = params["initial", "theta"]
    phi   = params["initial", "phi"]
    n1, n2, n3 = cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)
    xi = n1*X + n2*Y + n3*Z
    qexact = empty((mx, my, mz, 8))

    un = zeros((mx, my, mz))
    ut = 0.1 * sin(2.0*pi * (xi+t))
    ur = 0.1 * cos(2.0*pi * (xi+t))
    Bn = ones((mx, my, mz))
    Bt = 0.1 * sin(2.0*pi * (xi+t))
    Br = 0.1 * cos(2.0*pi * (xi+t))
    n1 = cos(phi) * cos(theta)
    n2 = sin(phi) * cos(theta)
    n3 = sin(theta)
    t1 = -sin(phi)
    t2 = cos(phi)
    t3 = 0.0
    r1 = -cos(phi)*sin(theta)
    r2 = -sin(phi)*sin(theta)
    r3 = cos(theta)

    qexact[:, :, :, 1 - 1] = 1.0
    qexact[:, :, :, 2 - 1] = n1*un + t1*ut + r1*ur
    qexact[:, :, :, 3 - 1] = n2*un + t2*ut + r2*ur
    qexact[:, :, :, 4 - 1] = n3*un + t3*ut + r3*ur
    qexact[:, :, :, 6 - 1] = n1*Bn + t1*Bt + r1*Br
    qexact[:, :, :, 7 - 1] = n2*Bn + t2*Bt + r2*Br
    qexact[:, :, :, 8 - 1] = n3*Bn + t3*Bt + r3*Br
    pexact = 0.1
    # We don't care about energy, because we don't test the accuracy
    # of that.
    # qexact[:, :, :, 5 - 1]
    #L1_error = sum(abs(q[:,:,:, 2 - 1]/q[:,:,:, 1 - 1] - qexact[:,:,:, 2 - 1]/qexact[:,:,:, 1 - 1])) / sum(abs(qexact[:,:,:, 2 - 1]/qexact[:,:,:, 1 - 1]))   #u1
    #L1_error = sum(abs(q[:,:,:, 4 - 1]/q[:,:,:, 1 - 1] - qexact[:,:,:, 4 - 1]/qexact[:,:,:, 1 - 1])) / sum(abs(qexact[:,:,:, 4 - 1]/qexact[:,:,:, 1 - 1]))   #u3
    #L1_error = sum(abs(q[:,:,:, 6 - 1] - qexact[:,:,:, 6 - 1])) / sum(abs(qexact[:,:,:, 6 - 1]))  #B1
    #L1_error = sum(abs(q[:,:,:, 8 - 1] - qexact[:,:,:, 8 - 1])) / sum(abs(qexact[:,:,:, 8 - 1]))  #B3
    #error_list.append(L1_error)
    #error = abs(q[:,:,:,0]-qexact[:,:,:,0])
    #error = abs(q[:,:,:, 6 - 1] - qexact[:,:,:, 6 - 1])  #B1
    error = abs(q[:,:,:, 1 - 1] - qexact[:,:,:, 1 - 1])   #rho
    error_list.append(max(error))

order_list = log2_adjacent_ratio(error_list)

print "Linfinity, error and order. rho"
print error_list
print order_list
