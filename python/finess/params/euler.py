"""Export parameter_list, accessor_list, check_list for Euler
equations.

Currently, only euler.gamma.
"""

from __future__ import absolute_import

def _parameters_accessors_checks():
    from finess.params import Parameter, CheckGreaterThan, Accessor
    parameters = []
    checks = []
    gamma = Parameter(variable_name = "gamma",
                      section = "euler",
                      name = "gamma",
                      type_ = "double",
                      default_value = 1.4
                      )
    parameters.append(gamma)
    checks.append(CheckGreaterThan(gamma, 0.0))

    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list =  \
    _parameters_accessors_checks()


