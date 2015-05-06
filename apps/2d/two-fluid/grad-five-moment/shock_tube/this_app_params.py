#section [initial]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
    """Routine to set initial conditions for a Riemann problem.
    
    TODO - this was written for MHD, this needs to be redone for 2D two-fluid!
    -DS
    """
     
    parameters  = []
    checks      = []
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
    
    u3l = Parameter(variable_name = "u3l",
                    section = "initial",
                    name = "u3l",
                    type_ = "double")
    parameters.append(u3l)
    
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
    
    B3l = Parameter(variable_name = "B3l",
                    section = "initial",
                    name = "B3l",
                    type_ = "double")
    parameters.append(B3l)
    
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
    
    u3r = Parameter(variable_name = "u3r",
                    section = "initial",
                    name = "u3r",
                    type_ = "double")
    parameters.append(u3r)
    
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
    
    B3r = Parameter(variable_name = "B3r",
                    section = "initial",
                    name = "B3r",
                    type_ = "double")
    parameters.append(B3r)
    
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

