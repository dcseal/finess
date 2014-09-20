from __future__ import absolute_import

def _parameters_accessors_checks():    
    from finess.params import Parameter, Accessor, Check, \
                              CheckGreaterEqual, CheckGreaterThan, \
    			  CheckOneOf, EnumParameterType
    
    parameters = []
    checks = []
    global_alpha = Parameter(variable_name = "global_alpha",
                         section = "finess",
                         name = "global_alpha",
                         type_ = "bool",
                         default_value = False)
    parameters.append(global_alpha)
    
    output_dir = Parameter(variable_name = "output_dir",
                           section = "finess",
    		       name = "output_dir",
    		       type_ = "std::string",
    		       default_value = "output")
    parameters.append(output_dir)		       

    ndims = Parameter(variable_name = "ndims",
                      section = "finess",
                      name = "ndims",
                      type_ = "int")
    parameters.append(ndims)
    checks.append(CheckOneOf(ndims, [1, 2, 3]))
    
    nout = Parameter(variable_name = "nout",
                     section = "finess",
                     name = "nout",
                     type_ = "int",
                     default_value = 1)
    parameters.append(nout)
    checks.append(CheckGreaterThan(nout, 0))
    
    tfinal = Parameter(variable_name = "tfinal",
                       section = "finess",
                       name = "tfinal",
                       type_ = "double")
    parameters.append(tfinal)
    checks.append(CheckGreaterEqual(tfinal, 0.0))
    
    initial_dt = Parameter(variable_name = "initial_dt",
                           section = "finess",
                           name = "initial_dt",
                           type_ = "double")
    parameters.append(initial_dt)
    checks.append(CheckGreaterThan(initial_dt, 0.0))
    
    max_dt = Parameter(variable_name = "max_dt",
                       section = "finess",
                       name = "max_dt",
                       type_ = "double")
    parameters.append(max_dt)
    checks.append(Check(cpp_code = """if(!(this->max_dt >= this->initial_dt))
        terminate("finess.max_dt should > finess.initial_dt");
    """))
    
    desired_cfl = Parameter(variable_name = "desired_cfl",
                            section = "finess",
                            name = "desired_cfl",
                            type_ = "double")
    parameters.append(desired_cfl)
    checks.append(CheckGreaterThan(desired_cfl, 0.0))
    
    max_cfl = Parameter(variable_name = "max_cfl",
                        section = "finess",
                        name = "max_cfl",
                        type_ = "double")
    parameters.append(max_cfl)
    checks.append(Check(cpp_code = """if(!(this->max_cfl >= this->desired_cfl))
        terminate("finess.max_cfl should > this->desired_cfl");
    """))
    
    #max number of time steps per call to FinSolvexxx
    nv = Parameter(variable_name = "nv",
                   section = "finess",
                   name = "nv",
                   type_ = "int")
    parameters.append(nv)
    checks.append(CheckGreaterThan(nv, 0))
    
    
    #TIGHTLY COUPLED with other C++ code
    #Need to change the other C++ code (RunFinpack)
    time_stepping_method = Parameter(variable_name = "time_stepping_method",
                                     section = "finess",
                                     name = "time_stepping_method",
                                     type_ = EnumParameterType(enum_scope_name = "TimeSteppingMethod",
                                                 string_enumerator_dict = \
                                                     {"Runge-Kutta": "RK",
                                                      "Lax-Wendroff": "LxW",
                                                      "Multiderivative": "MD",
                                                      "User-Defined": "USER_DEFINED"} ))
    parameters.append(time_stepping_method)
    
    space_order = Parameter(variable_name = "space_order",
                            section = "finess",
                            name = "space_order",
                            type_ = "int")
    parameters.append(space_order)
    checks.append(CheckOneOf(space_order, [1, 3, 5, 7, 9, 11]))
    
    time_order = Parameter(variable_name = "time_order",
                           section = "finess",
                           name = "time_order",
                           type_ = "int")
    parameters.append(time_order)
    checks.append(CheckOneOf(time_order, [1, 2, 3, 4, 5]))
    
    #use_limiter removed
    
    verbosity = Parameter(variable_name = "verbosity",
                          section = "finess",
                          name = "verbosity",
                          type_ = "int")
    parameters.append(verbosity)
    checks.append(CheckOneOf(verbosity, [0, 1]))
    
    mcapa = Parameter(variable_name = "mcapa",
                      section = "finess",
                      name = "mcapa",
                      type_ = "int")
    parameters.append(mcapa)
    checks.append(CheckGreaterEqual(mcapa, 0))
    
    maux = Parameter(variable_name = "maux",
                     section = "finess",
                     name = "maux",
                     type_ = "int")
    parameters.append(maux)
    checks.append(Check(cpp_code="""if(!(this->maux >= this->mcapa))
        terminate("finess.maux should >= finess.mcapa");
    """))
    
    source_term = Parameter(variable_name = "source_term",
                            section = "finess",
                            name = "source_term",
                            type_ = "bool")
    parameters.append(source_term)
    
    meqn = Parameter(variable_name = "meqn",
                     section = "finess",
                     name = "meqn",
                     type_ = "int")
    parameters.append(meqn)
    checks.append(CheckGreaterEqual(meqn, 1))

    return parameters, map(Accessor, parameters), checks
    
parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

