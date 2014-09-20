from __future__ import absolute_import

def _parameters_accessors_checks():
    from finess.params import Parameter, DerivedParameter, Accessor,\
         Check, CheckGreaterEqual, CheckGreaterThan
    parameters = []
    checks = []
    mx = Parameter(variable_name = "mx",
                   section = "grid",
                   name = "mx",
                   type_ = "int")
    parameters.append(mx)
    checks.append(CheckGreaterThan(mx, 0))
    
    my = Parameter(variable_name = "my",
                   section = "grid",
                   name = "my",
                   type_ = "int")
    parameters.append(my)
    checks.append(CheckGreaterThan(my, 0))
    mbc = Parameter(variable_name = "mbc",
                    section = "grid",
                    name = "mbc",
                    type_ = "int")
    parameters.append(mbc)
    checks.append(CheckGreaterEqual(mbc, 0))

    xlow = Parameter(variable_name = "xlow",
                     section = "grid",
		 name = "xlow",
		 type_ = "double")
    parameters.append(xlow)
    
    xhigh = Parameter(variable_name = "xhigh",
                      section = "grid",
		  name = "xhigh",
		  type_ = "double")
    parameters.append(xhigh)
    checks.append(Check(cpp_code = """if(!(xhigh > xlow))
        terminate("grid.xhigh should > grid.xlow.");"""))
    
    
    ylow = Parameter(variable_name = "ylow",
                     section = "grid",
		 name = "ylow",
		 type_ = "double")
    parameters.append(ylow)
    
    yhigh = Parameter(variable_name = "yhigh",
                      section = "grid",
		  name = "yhigh",
		  type_ = "double")
    parameters.append(yhigh)
    checks.append(Check(cpp_code = """if(!(yhigh > ylow))
        terminate("grid.yhigh should > grid.ylow.");"""))
    
    dx = DerivedParameter(variable_name = "dx",
                          type_ = "double",
		      defining_expression_in_cpp = """(this->xhigh - this->xlow) / this->mx""")
    parameters.append(dx)
    
    dy = DerivedParameter(variable_name = "dy",
                          type_ = "double",
		      defining_expression_in_cpp = """(this->yhigh - this->ylow) / this->my""")
    parameters.append(dy)
    
    return parameters, map(Accessor, parameters), checks



parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

