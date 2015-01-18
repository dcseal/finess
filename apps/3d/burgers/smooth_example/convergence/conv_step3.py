# Generate exact solutions for Burgers' equation


def generate_exact_solution(parameters_ini_filename,
                           q0,  max_abs_sum_partial_deriv_q0):
    from finess.params.util import read_params
    
    from generate_iniparams import parameter_list
    from burger_exact import burger_exact_3d
    
    params = read_params(parameters_ini_filename, parameter_list)
    output_dir = params['finess', 'output_dir']
    import os
    if not os.path.isdir(output_dir):
        # If output_dir exists and is not a dir, then this fails,
        # which is the correct behavior
        os.mkdir(output_dir)
    output_filename = os.path.join(output_dir, 'qexact.silo')
    
    mx    = params['grid', 'mx']
    my    = params['grid', 'my']
    mz    = params['grid', 'mz']
    xlow  = params['grid', 'xlow']
    xhigh = params['grid', 'xhigh']
    ylow  = params['grid', 'ylow']
    yhigh = params['grid', 'yhigh']
    zlow  = params['grid', 'zlow']
    zhigh = params['grid', 'zhigh']
    
    dx = (xhigh-xlow) / float(mx)
    dy = (yhigh-ylow) / float(my)
    dz = (zhigh-zlow) / float(mz)
    tfinal = params['finess', 'tfinal']
    
    from numpy import empty
    exact_solution = empty((mx, my, mz))
    
    for i in range(mx):
        for j in range(my):
            for k in range(mz):
                x = xlow + (i+0.5)*dx
                y = ylow + (j+0.5)*dy
                z = zlow + (k+0.5)*dz
                q = burger_exact_3d(x, y, z, tfinal, q0, max_abs_sum_partial_deriv_q0)
                exact_solution[i, j, k] = q
    
    from pyvisfile import silo
    from numpy import linspace, asarray
    f = silo.SiloFile(output_filename, mode=silo.DB_CLOBBER)
    coord = [linspace(xlow + 0.5*dx, xhigh - 0.5*dx, mx),
             linspace(ylow + 0.5*dy, yhigh - 0.5*dy, my), 
             linspace(zlow + 0.5*dz, zhigh - 0.5*dz, mz)]
    options = dict()
    options[silo.DBOPT_DTIME] = tfinal
    f.put_quadmesh("quadmesh", coord, optlist=options)
    
    f.put_quadvar1("q", "quadmesh",
                   asarray(exact_solution, order="F"),
                   exact_solution.shape,
                   centering=silo.DB_NODECENT)
    f.close()


parameters_ini_filename_list = \
    ["parameters%02d.ini" % i for i in [0, 1, 2, 3, 4]]

from math import sin, pi
q0 = lambda x, y, z: 0.25 + sin(pi*x) * sin(pi*y) * sin(pi*z)
max_abs_sum_partial_deriv_q0 = 3.0 * pi

for parameters_ini_filename in parameters_ini_filename_list:
    print "Now processing %s." % parameters_ini_filename
    generate_exact_solution(parameters_ini_filename,
                            q0, max_abs_sum_partial_deriv_q0)

                
