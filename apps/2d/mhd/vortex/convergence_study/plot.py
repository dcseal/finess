# -*- coding: utf-8 -*-
"""
Created on Sun Dec  7 18:27:05 2014

@author: xiaofeng
"""


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

output_dir = "output04"
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
gamma = params["mhd", "gamma"]

tfinal, q, aux = finess.dim2.read_qa(params, nout)

from numpy import sqrt

rho = q[:, :, 1 - 1]
u1 = q[:, :, 2 - 1] / rho
u2 = q[:, :, 3 - 1] / rho
u3 = q[:, :, 4 - 1] / rho
E = q[:, :, 5 - 1]
B1 = q[:, :, 6 - 1]
B2 = q[:, :, 7 - 1]
B3 = q[:, :, 8 - 1]

upertubnorm = sqrt((u1-1)**2 + (u2-1)**2 + u3**2)
Bnorm = sqrt(B1**2 + B2**2 + B3**2)

p = (gamma - 1) * (E - 0.5 * rho * (u1**2 + u2**2 + u3**2) - 0.5 * (B1**2 + B2**2 + B3**2))

import finess.viz.dim2

X, Y = finess.viz.dim2.meshgrid(params)

import pylab
#pylab.plot(X, Y, logrho)

nContourLevels = 30

pylab.figure(1)
pylab.title('$|\mathbf{u} - (1,1)|$, $t=%(tfinal)f$' % {"tfinal": tfinal})
pylab.contour(X, Y, upertubnorm, nContourLevels)
pylab.axes().set_aspect('equal')

pylab.figure(2)
pylab.title('$|\mathbf{B}|$, $t=%(tfinal)f$' % {"tfinal": tfinal})
pylab.contour(X, Y, Bnorm, nContourLevels)
pylab.axes().set_aspect('equal')

pylab.figure(3)
pylab.title('$p$, $t=%(tfinal)f$' % {"tfinal": tfinal})
pylab.contour(X, Y, p, nContourLevels)
pylab.axes().set_aspect('equal')

pylab.figure(4)
pylab.title(r"""$\rho$, $t=%(tfinal)f$""" % {"tfinal": tfinal})
pylab.contour(X, Y, rho, nContourLevels)
pylab.axes().set_aspect('equal')


pylab.show()

pylab.gradient

