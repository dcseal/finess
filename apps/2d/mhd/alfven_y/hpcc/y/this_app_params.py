#section [initial]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
     
    parameters = []
    checks = []
    
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

