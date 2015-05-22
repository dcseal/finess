def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
     
    parameters = []
    checks = []

    keep_pressure = Parameter(variable_name = "keep_pressure",
                      section = "mhd",
                      name = "keep_pressure",
                      type_ = "bool")
    parameters.append(keep_pressure)

    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

