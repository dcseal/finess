#section [initial]
def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
     
    parameters  = []
    checks      = []

    x0 = Parameter(variable_name = "x0",
                      section = "euler",
                      name = "x0",
                      type_ = "double",
                      default_value = 0.0
                      )
    parameters.append(x0)
#   checks.append(CheckGreaterThan(x0, 0.0))

    riemann_problem_no = Parameter(variable_name = "riemann_problem_number",
                    section = "euler",
                    name = "riemann_problem_number",
                    type_ = "int")
    parameters.append(riemann_problem_no)
    checks.append(CheckOneOf(riemann_problem_no, [0, 1, 2, 3, 4, 5]))
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

