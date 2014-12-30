#section [initialparams]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
    """Set user-defined parameters for this application.

    This routine saves values used for a Riemann problem.  They are saved in
    rho[l,r], u1[l,r], u2[l,r], u3[l,r] and p[l,r] for left, and right states.
    """

     
    parameters = []
    checks = []
    rhol = Parameter(variable_name = "rhol",
                     section = "initialparams",
                     name = "rhol",
                     type_ = "double")
    parameters.append(rhol)
    checks.append(CheckGreaterThan(rhol, 0.0))
    
    unl = Parameter(variable_name = "u1l",
                    section = "initialparams",
                    name = "u1l",
                    type_ = "double")
    parameters.append(unl)
    
    utl = Parameter(variable_name = "u2l",
                    section = "initialparams",
                    name = "u2l",
                    type_ = "double")
    parameters.append(utl)
    
    u3l = Parameter(variable_name = "u3l",
                    section = "initialparams",
                    name = "u3l",
                    type_ = "double")
    parameters.append(u3l)
    
    pl = Parameter(variable_name = "pl",
                   section = "initialparams",
                   name = "pl",
                   type_ = "double")
    parameters.append(pl)
    checks.append(CheckGreaterThan(pl, 0.0))
    
    rhor = Parameter(variable_name = "rhor",
                     section = "initialparams",
                     name = "rhor",
                     type_ = "double")
    parameters.append(rhor)
    checks.append(CheckGreaterThan(rhor, 0.0))
    
    unr = Parameter(variable_name = "u1r",
                    section = "initialparams",
                    name = "u1r",
                    type_ = "double")
    parameters.append(unr)
    
    utr = Parameter(variable_name = "u2r",
                    section = "initialparams",
                    name = "u2r",
                    type_ = "double")
    parameters.append(utr)
    
    u3r = Parameter(variable_name = "u3r",
                    section = "initialparams",
                    name = "u3r",
                    type_ = "double")
    parameters.append(u3r)
    
    pr = Parameter(variable_name = "pr",
                   section = "initialparams",
                   name = "pr",
                   type_ = "double")
    parameters.append(pr)
    checks.append(CheckGreaterThan(pr, 0.0))
    
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

