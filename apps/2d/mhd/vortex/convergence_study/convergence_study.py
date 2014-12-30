# coding: utf-8



def get_dict_output_dir_to_parameters_ini_dump_filename():
    import os
    dir_ = '.'
    output_dir_list = sorted([output_dir for output_dir in os.listdir(dir_) if output_dir.startswith('output')])
    ret = {}
    for output_dir in output_dir_list:
        with open(os.path.join(output_dir, 'parameters_ini_filename')) as f:
            parameters_ini_filename = list(f)[0].rstrip()
        ret[output_dir] = parameters_ini_filename + '.dump'
    return ret

        
dict_output_dir_to_parameters_ini_dump = get_dict_output_dir_to_parameters_ini_dump_filename()


import finess.util
import finess.params.util
import finess.dim2
import generate_iniparams


#     q(:, :, i - 1):
#     * i = 1: mass
#     * i = 2: momentum-1
#     * i = 3: momentum-2
#     * i = 4: momentum-3
#     * i = 5: energy
#     * i = 6: B1
#     * i = 7: B2
#     * i = 8: B3

import finess.viz.dim2

def L1_error_list(output_dir_list):
    global debug_rho
    
    from numpy import sum, abs, sqrt, pi, exp, empty_like, max
    import finess.viz.dim2
    error_list = []
    for output_dir in output_dir_list:
        parameters_ini_dump_filename = dict_output_dir_to_parameters_ini_dump[output_dir]
        import os.path
        params = finess.params.util.read_params(os.path.join(output_dir, parameters_ini_dump_filename), generate_iniparams.parameter_list)
        xlow = params['grid', 'xlow']
        xhigh = params['grid', 'xhigh']
        ylow = params['grid', 'ylow']
        yhigh = params['grid', 'yhigh']
        mx = params['grid', 'mx']
        my = params['grid', 'my']
        dx = (xhigh - xlow) / float(mx)
        dy = (yhigh - ylow) / float(my)
        nout = params['finess', 'nout']
        gamma = params['mhd', 'gamma']
        kappa = params['vortex', 'kappa']
        mu = params['vortex', 'mu']

        tfinal, q, aux = finess.dim2.read_qa(params, nout)
        print "tfinal: ", tfinal
        X, Y = finess.viz.dim2.meshgrid(params)
        
        
        rho = q[:, :, 1 - 1]
        debug_rho = rho
        u1 = q[:, :, 2 - 1] / q[:, :, 1 - 1]
        u2 = q[:, :, 3 - 1] / q[:, :, 1 - 1]
        u3 = q[:, :, 4 - 1] / q[:, :, 1 - 1]
        E = q[:, :, 5 - 1]
        B1 = q[:, :, 6 - 1]
        B2 = q[:, :, 7 - 1]
        B3 = q[:, :, 8 - 1]
        p = (gamma - 1) * (E - 0.5*rho*(u1**2 + u2**2 + u3**2) \
                             - 0.5*(B1**2 + B2**2 + B3**2))
        
        x_displacement = u1_mean * tfinal
        y_displacement = u2_mean * tfinal
        # 'before displaced'
        Xbd = X - x_displacement
        Ybd = Y - y_displacement
        r2bd = Xbd**2 + Ybd**2
        rho_exact = empty_like(rho)
        rho_exact.fill(1.0)
        u1_exact= u1_mean + kappa/(2.0*pi) * exp(0.5*(1.0-r2bd)) * (-Ybd)
        u2_exact= u2_mean + kappa/(2.0*pi) * exp(0.5*(1.0-r2bd)) * Xbd
        u3_exact = 0.0
        B1_exact = B1_mean + mu/(2.0*pi) * exp(0.5*(1.0-r2bd)) * (-Ybd)
        B2_exact = B2_mean + mu/(2.0*pi) * exp(0.5*(1.0-r2bd)) * Xbd
        B3_exact = 0.0
        p_exact = p_mean + 1.0/(8.0*pi**2) * exp(1-r2bd) \
                          *(mu**2 * (1-r2bd) - kappa**2) 

        L1_error_rho = sum(abs(rho - rho_exact))
        #L1_error_unorm = sum(abs(sqrt(u1**2 + u2**2 + u3**2) - sqrt(u1_exact**2 + u2_exact**2 + u3_exact**2)))
        #L1_error_B = sum(abs(sqrt((B1-B1_exact)**2 + (B2-B2_exact)**2 + (B3-B3_exact)**2)))
        L1_error_p = sum(abs(p - p_exact))

        L1_rho_exact = sum(abs(rho_exact))
        #L1_unorm_exact = sum(sqrt(u1_exact**2 + u2_exact**2 + u3_exact**2))
        #L1_Bnorm_exact = sum(sqrt(B1_exact**2 + B2_exact**2 + B3_exact**2))
        L1_p_exact = sum(abs(p_exact))
        
        L1_error_B1 = sum(abs(B1-B1_exact))
        L1_B1_exact = sum(abs(B1_exact))
        L1_error_B2 = sum(abs(B2-B2_exact))
        L1_B2_exact = sum(abs(B2_exact))
        
        L1_error_u1 = sum(abs(u1-u1_exact))
        L1_u1_exact = sum(abs(u1_exact))
        L1_error_u2 = sum(abs(u2-u2_exact))
        L1_u2_exact = sum(abs(u2_exact))        

        
        Linfinity_error_rho = max(abs(rho - rho_exact))
        Linfinity_error_p = max(abs(p - p_exact))
        
        Linfinity_error_B1 = max(abs(B1-B1_exact))        
        Linfinity_error_B2 = max(abs(B2-B2_exact))        
        
        Linfinity_error_u1 = max(abs(u1-u1_exact))        
        Linfinity_error_u2 = max(abs(u2-u2_exact))            
        
        #delta = Linfinity_error_u1
        #delta = L1_error_rho/L1_rho_exact  #good
        #delta = L1_error_B1/L1_B1_exact
        #delta = L1_error_B2/L1_B2_exact
        #delta = L1_error_u2/L1_u2_exact
        #delta = L1_error_p/L1_p_exact
        delta = Linfinity_error_rho
        #delta = Linfinity_error_p
        error_list.append(delta)
    return error_list


def log2_adjacent_ratio(error_list):
    order_list = []
    from numpy import log2
    for i in range(len(error_list) - 1):
        order_list.append(log2(error_list[i] / error_list[i+1]))
    return order_list

#def L1_A_error_list(output_dir_list):
#    from numpy import exp, pi    
#    import finess.viz.dim2
#    error_list = []
#    for output_dir in output_dir_list:
#        parameters_ini_dump_filename = dict_output_dir_to_parameters_ini_dump[output_dir]
#        import os.path
#        params = finess.params.util.read_params(os.path.join(output_dir, parameters_ini_dump_filename), generate_iniparams.parameter_list)
#        xlow = params['grid', 'xlow']
#        xhigh = params['grid', 'xhigh']
#        ylow = params['grid', 'ylow']
#        yhigh = params['grid', 'yhigh']
#        mx = params['grid', 'mx']
#        my = params['grid', 'my']
#        dx = (xhigh - xlow) / float(mx)
#        dy = (yhigh - ylow) / float(my)
#        nout = params['finess', 'nout']
#        mu = params['vortex', 'mu']
#        tfinal, q, aux = finess.dim2.read_qa(params, nout)
#        print "tfinal: ", tfinal
#        X, Y = finess.viz.dim2.meshgrid(params)
#        
#        
#        x_displacement = u1_mean * tfinal
#        y_displacement = u2_mean * tfinal
#        # 'before displaced'
#        Xbd = X - x_displacement
#        Ybd = Y - y_displacement
#        r2bd = Xbd**2 + Ybd**2
#
#        A_exact = mu/(2.0*pi) * exp(0.5*(1.0-r2bd))
#        A = aux[:, :, 1 - 1]
#        
#        from numpy import sum, abs
#        
#        L1_A_exact = sum(abs(A_exact))
#        L1_A_error = sum(abs(A - A_exact))
#        print 'L1 A_exact:', L1_A_exact
#        print 'L1 A error:', L1_A_error
#        delta = L1_A_error / L1_A_exact
#        print 'delta:', delta
#        error_list.append(delta)
#    return error_list



u1_mean = 1.0
u2_mean = 1.0
B1_mean = 0.0
B2_mean = 0.0
p_mean = 1.0

output_dir_list = ['output%(i)02d' % {'i': i} for i in [0, 1, 2, 3, 4]]
error_list = L1_error_list(output_dir_list)
order_list = log2_adjacent_ratio(error_list)
print order_list
print error_list

#
#
#A_error_list = L1_A_error_list(output_dir_list)
#A_order_list = log2_adjacent_ratio(A_error_list)
#print 'A:'
#print A_order_list
#print A_error_list
#
