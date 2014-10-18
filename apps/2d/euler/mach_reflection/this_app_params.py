#section [initial]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
     
    parameters = []
    checks = []

    parameters.append(Parameter(variable_name = "x0",
                      section = "euler",
                      name = "x0",
                      type_ = "double"))
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

