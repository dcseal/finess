#section [initial]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
     
    parameters = []
    checks = []
    rhol = Parameter(variable_name = "rhol",
                     section = "initial",
                     name = "rhol",
                     type_ = "double")
    parameters.append(rhol)
    checks.append(CheckGreaterThan(rhol, 0.0))
    
    unl = Parameter(variable_name = "unl",
                    section = "initial",
                    name = "unl",
                    type_ = "double")
    parameters.append(unl)
    
    utl = Parameter(variable_name = "utl",
                    section = "initial",
                    name = "utl",
                    type_ = "double")
    parameters.append(utl)
    
    url = Parameter(variable_name = "url",
                    section = "initial",
                    name = "url",
                    type_ = "double")
    parameters.append(url)
    
    pl = Parameter(variable_name = "pl",
                   section = "initial",
                   name = "pl",
                   type_ = "double")
    parameters.append(pl)
    checks.append(CheckGreaterThan(pl, 0.0))
    
    Bnl = Parameter(variable_name = "Bnl",
                    section = "initial",
                    name = "Bnl",
                    type_ = "double")
    parameters.append(Bnl)
    
    Btl = Parameter(variable_name = "Btl",
                    section = "initial",
                    name = "Btl",
                    type_ = "double")
    parameters.append(Btl)
    
    Brl = Parameter(variable_name = "Brl",
                    section = "initial",
                    name = "Brl",
                    type_ = "double")
    parameters.append(Brl)
    
    rhor = Parameter(variable_name = "rhor",
                     section = "initial",
                     name = "rhor",
                     type_ = "double")
    parameters.append(rhor)
    checks.append(CheckGreaterThan(rhor, 0.0))
    
    unr = Parameter(variable_name = "unr",
                    section = "initial",
                    name = "unr",
                    type_ = "double")
    parameters.append(unr)
    
    utr = Parameter(variable_name = "utr",
                    section = "initial",
                    name = "utr",
                    type_ = "double")
    parameters.append(utr)
    
    urr = Parameter(variable_name = "urr",
                    section = "initial",
                    name = "urr",
                    type_ = "double")
    parameters.append(urr)
    
    pr = Parameter(variable_name = "pr",
                   section = "initial",
                   name = "pr",
                   type_ = "double")
    parameters.append(pr)
    checks.append(CheckGreaterThan(pr, 0.0))
    
    Bnr = Parameter(variable_name = "Bnr",
                    section = "initial",
                    name = "Bnr",
                    type_ = "double")
    parameters.append(Bnr)
    
    Btr = Parameter(variable_name = "Btr",
                    section = "initial",
                    name = "Btr",
                    type_ = "double")
    parameters.append(Btr)
    
    Brr = Parameter(variable_name = "Brr",
                    section = "initial",
                    name = "Brr",
                    type_ = "double")
    parameters.append(Brr)
    
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

