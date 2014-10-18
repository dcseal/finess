"""Export parameter_list, accessor_list, check_list for WENO
reconstructions.


"""

from __future__ import absolute_import

def _parameters_accessors_checks():
    from finess.params import Parameter, EnumParameterType, \
         CheckGreaterEqual, CheckGreaterThan, Accessor
    parameters = []
    checks = []

    weno_version = Parameter(variable_name = "weno_version",
                       section = "weno",
    		   name = "weno_version",
    		   type_ = EnumParameterType( \
    		             enum_scope_name = "WenoVersion",
    			     string_enumerator_dict = \
    			         {"FD": "FD",
    				  "JS": "JS",
    				  "Z": "Z"}),
    	           default_value = "JS")
    parameters.append(weno_version)
    
    power_param = Parameter(variable_name = "power_param",
                           section = "weno",
    		       name = "power_param",
    		       type_ = "double",
    		       default_value = 2.0)
    parameters.append(power_param)
    
    
    
    alpha_scaling = Parameter(variable_name = "alpha_scaling",
                              section = "weno",
                              name = "alpha_scaling",
                              type_ = "double",
                              default_value = 1.1)
    parameters.append(alpha_scaling)
    checks.append(CheckGreaterEqual(alpha_scaling, 1.0))
    
    
    epsilon = Parameter(variable_name = "epsilon",
                        section = "weno",
                        name = "epsilon",
                        type_ = "double",
                        default_value = 1e-6)
    parameters.append(epsilon)
    checks.append(CheckGreaterThan(epsilon, 0.0))
    
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()
