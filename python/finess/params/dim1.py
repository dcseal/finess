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
    
    
    dx = DerivedParameter(variable_name = "dx",
                          type_ = "double",
		      defining_expression_in_cpp = """(this->xhigh - this->xlow) / this->mx""")
    parameters.append(dx)
    
    return parameters, map(Accessor, parameters), checks



parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()


def starter_pac():
    from finess.params import finess, weno, dim1, append_pac_from_module
    
    pac = ([], [], [])
    
    append_pac_from_module(pac, finess)
    append_pac_from_module(pac, weno)   
    append_pac_from_module(pac, dim1)
       
    return pac 

