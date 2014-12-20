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

output_dir = "output"
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

tfinal, q, aux = finess.dim2.read_qa(params, 5)

from numpy import log, gradient, sqrt
logrho = log(q[:, :, 1 - 1])
gradX, gradY = gradient(logrho, dx, dy)
normgradlogrho = sqrt(gradX**2 + gradY**2)

import finess.viz.dim2

X, Y = finess.viz.dim2.meshgrid(params)

import pylab
#pylab.plot(X, Y, logrho)

pylab.pcolormesh(X, Y, normgradlogrho, cmap = "gist_stern")

pylab.show()

pylab.gradient

