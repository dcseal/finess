"""Export parameter_list, accessor_list, check_list for MHD 
equations.

Currently, only mhd.gamma.
"""

from __future__ import absolute_import

def _parameters_accessors_checks():
    from finess.params import Parameter, CheckGreaterThan, Accessor
    parameters = []
    checks = []
    gamma = Parameter(variable_name = "gamma",
                      section = "mhd",
                      name = "gamma",
                      type_ = "double")
    parameters.append(gamma)
    checks.append(CheckGreaterThan(gamma, 0.0))
    
    constrained_transport = Parameter(variable_name = "constrained_transport",
                                      section = "mhd",
                                      name = "constrained_transport",
                                      type_ = "bool",
                                      default_value = "true")
    parameters.append(constrained_transport)

    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list =  \
    _parameters_accessors_checks()


