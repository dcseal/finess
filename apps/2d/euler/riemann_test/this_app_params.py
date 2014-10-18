#section [initial]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
     
    parameters = []
    checks = []

    OPT = Parameter(variable_name = "OPT",
                    section = "euler",
                    name = "OPT",
                    type_ = "int")
    parameters.append(OPT)
    checks.append(CheckOneOf(OPT, [1, 2]))
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

