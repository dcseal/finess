# Generate exact solutions for Burgers' equation


def log2_adjacent_ratio(error_list):
    order_list = []
    from numpy import log2
    for i in range(len(error_list) - 1):
        order_list.append(log2(error_list[i] / error_list[i+1]))
    return order_list


parameters_ini_filename_list = \
    ["output%(i)02d/parameters%(i)02d.ini.dump" % {"i": i} for i in [0, 1, 2]]


error_list = []
error_array_list = []

B1_list = []
B2_list = []
B3_list = []

from numpy import abs, max, sum, zeros, ones, argmax, unravel_index
from numpy import sin, cos, pi
for parameters_ini_filename in parameters_ini_filename_list:
    from numpy import empty
    from finess.params.util import read_params
    from generate_iniparams import parameter_list
    from finess.dim3 import read_qa
    params = read_params(parameters_ini_filename, parameter_list)
    print "Processing %s" % parameters_ini_filename
    nout = params['finess', 'nout']
    t, q, aux = read_qa(params, nout)
    
    B1_list.append(q[:, :, :, 6-1])
    B2_list.append(q[:, :, :, 7-1])
    B3_list.append(q[:, :, :, 8-1])
    from finess.viz.dim3 import meshgrid
    mx = params["grid", "mx"]
    my = params["grid", "my"]
    mz = params["grid", "mz"]

    X, Y, Z = meshgrid(params)
    print X.shape, Y.shape, Z.shape
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
    t1 = -sin(phi)
    t2 = cos(phi)
    t3 = 0.0
    r1 = -cos(phi)*sin(theta)
    r2 = -sin(phi)*sin(theta)
    r3 = cos(theta)
    
    B1av = cos(phi) * cos(theta)
    B2av = sin(phi) * cos(theta)
    B3av = sin(theta)
    
    aexact = empty((mx, my, mz, 3))
    aexact[:, :, :, 1 - 1] = Z * B2av - 1.0/(20.0*pi) * sin(phi) * sin(2.0*pi*(xi+t))
    aexact[:, :, :, 2 - 1] = X * B3av + 1.0/(20.0*pi) * cos(phi) * sin(2.0*pi*(xi+t))
    aexact[:, :, :, 3 - 1] = Y * B1av + 1.0/(20.0*pi*cos(theta)) * cos(2.0*pi*(xi+t))
    
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
    #error = abs(q[:,:,:, 1 - 1] - qexact[:,:,:, 1 - 1])   #rho
    #error = abs(q[:,:,:, 7 - 1] - qexact[:,:,:, 7 - 1])  #B2
    #error = abs(q[:,:,:, 8 - 1] - qexact[:,:,:, 8 - 1])  #B3
    #error = abs(q[:, :, :, 2 - 1]/q[:, :, :, 1-1] - qexact[:,:,:, 2-1]/qexact[:,:,:,1-1])  #u1
    #error = abs(q[:, :, :, 3 - 1]/q[:, :, :, 1-1] - qexact[:,:,:, 3-1]/qexact[:,:,:,1-1])  #u2
    error = abs(q[:, :, :, 4 - 1]/q[:, :, :, 1-1] - qexact[:,:,:, 4-1]/qexact[:,:,:,1-1])  #u3
    #error = abs(aux[:, :, :, 1 - 1] - aexact[:, :, :, 1 - 1])  #A1
    #error = abs(aux[:, :, :, 2 - 1] - aexact[:, :, :, 2 - 1])  #A2
    #error = abs(aux[:, :, :, 3 - 1] - aexact[:, :, :, 3 - 1])  #A3
    flattened_index_max_error = argmax(error)
    index_max_error = unravel_index(flattened_index_max_error, error.shape)
    print error.shape, index_max_error
    error_array_list.append(error)
    error_list.append(max(error))

order_list = log2_adjacent_ratio(error_list)

print "Linfinity, error and order. u3"
print error_list
print order_list
