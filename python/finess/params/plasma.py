"""Export parameter_list, accessor_list, check_list for two-fluid plasma equations.
"""

from __future__ import absolute_import

def _parameters_accessors_checks():

    from finess.params import Parameter, DerivedParameter, Accessor,\
         Check, CheckGreaterEqual, CheckGreaterThan

    parameters = []
    checks     = []

    # Name of plasma model (for now, only g05 is implemented)
    output_dir = Parameter(variable_name = "model_name",
                   section = "plasma",
    		       name    = "model_name",
    		       type_   = "std::string",
    		       default_value = "g05")
    parameters.append(output_dir)		       

    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list =  \
    _parameters_accessors_checks()


