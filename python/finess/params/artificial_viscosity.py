"""Export parameter_list, accessor_list, check_list for artificial
viscosity.

Reference: doi:10.1016/j.jcp.2014.03.001
"""

from __future__ import absolute_import

def _parameters_accessors_checks():
    from finess.params import Parameter, CheckGreaterEqual, Accessor
    parameters = []
    checks = []
    av_nu = Parameter(variable_name = "av_nu",
                      section = "av",
                      name = "nu",
                      type_ = "double")
    parameters.append(av_nu)
    checks.append(CheckGreaterEqual(av_nu, 0.0))

    av_delta = Parameter(variable_name = "av_delta",
                      section = "av",
                      name = "delta",
                      type_ = "double")
    parameters.append(av_delta)
    checks.append(CheckGreaterEqual(av_delta, 0.0))
    
    av_epsilon = Parameter(variable_name = "av_epsilon",
                      section = "av",
                      name = "epsilon",
                      type_ = "double")
    parameters.append(av_epsilon)
    checks.append(CheckGreaterEqual(av_epsilon, 0.0))
    
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list =  \
    _parameters_accessors_checks()


