#section [initial]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
     
    parameters = []
    checks = []
    
    checks.append(Check(cpp_code = """if(this->maux != 1)
        terminate("finess.maux should be 1 for 2D Constraint Transport.");"""))

    kappa = Parameter(variable_name = "kappa",
                      section = "vortex",
                      name = "kappa",
                      type_ = "double")
    parameters.append(kappa)

    mu = Parameter(variable_name = "mu",
                   section = "vortex",
                   name = "mu",
                   type_ = "double")
    parameters.append(mu)

    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

