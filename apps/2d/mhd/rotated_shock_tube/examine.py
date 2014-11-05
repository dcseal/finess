# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/xiaofeng/.spyder2/.temp.py
"""

import finess.dim2
import finess.params.util
import generate_iniparams

params = finess.params.util.read_params('output/parameters.ini.dump', 
                                        generate_iniparams.parameter_list)


t, q, aux = finess.dim2.read_qa(params, 60)
