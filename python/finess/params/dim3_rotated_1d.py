# Reference: arXiv:1007.2606, Section 6.2.1

from __future__ import absolute_import

def _parameters_accessors_checks():
    from finess.params import Parameter, DerivedParameter, Accessor,\
         Check, CheckGreaterEqual, CheckGreaterThan, CheckLessThan
    from math import pi

    parameters = []
    checks     = []

    # Grid information
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

    mz = Parameter(variable_name = "mz",
                   section = "grid",
                   name = "mz",
                   type_ = "int")
    parameters.append(mz)
    checks.append(CheckGreaterThan(mz, 0))

    # Number of ghost cells
    mbc = Parameter(variable_name = "mbc",
                    section = "grid",
                    name = "mbc",
                    type_ = "int")
    parameters.append(mbc)
    checks.append(CheckGreaterEqual(mbc, 0))

    # Domain information
    theta = Parameter(variable_name = 'theta',
                      section = 'initial',
                      name = 'theta',
                      type_ = 'double')
    parameters.append(theta)
    checks.append(CheckGreaterThan(theta, 0))
    checks.append(CheckLessThan(theta, pi / 2))

    phi = Parameter(variable_name = 'phi',
                      section = 'initial',
                      name = 'phi',
                      type_ = 'double')
    parameters.append(phi)
    checks.append(CheckGreaterThan(phi, 0))
    checks.append(CheckLessThan(phi, pi / 2))

   
    xlow = DerivedParameter(variable_name = "xlow",
                            type_ = "double",
                            defining_expression_in_cpp = "0.0",
                            dump_to_section = 'grid',
                            dump_to_name = 'xlow')
    parameters.append(xlow)
    
    xhigh = DerivedParameter(variable_name = "xhigh",
                             type_ = "double",
                             defining_expression_in_cpp = "1 / std::cos(this->theta) / std::cos(this->phi)",
                             dump_to_section = 'grid',
                             dump_to_name = 'xhigh')
    parameters.append(xhigh)
    
    ylow = DerivedParameter(variable_name = "ylow",
                            type_ = "double",
                            defining_expression_in_cpp = "0.0",
                            dump_to_section = 'grid',
                            dump_to_name = 'ylow')
    parameters.append(ylow)
   
    yhigh = DerivedParameter(variable_name = "yhigh",
                             type_ = "double",
                             defining_expression_in_cpp = "1 / std::cos(this->theta) / std::sin(this->phi)",
                             dump_to_section = 'grid',
                             dump_to_name = 'yhigh')
    parameters.append(yhigh)
   
    zlow = DerivedParameter(variable_name = "zlow",
                            type_ = "double",
                            defining_expression_in_cpp = "0.0",
                            dump_to_section = 'grid',
                            dump_to_name = 'zlow')
    parameters.append(zlow)
    
    zhigh = DerivedParameter(variable_name = "zhigh",
                             type_ = "double",
                             defining_expression_in_cpp = "1 / std::sin(this->theta)",
                             dump_to_section = 'grid',
                             dump_to_name = 'zhigh')
    parameters.append(zhigh)

    # Derived parameters (from grid)
    dx = DerivedParameter(variable_name = "dx",
                          type_ = "double",
		      defining_expression_in_cpp = """(this->xhigh - this->xlow) / this->mx""")
    parameters.append(dx)
    
    dy = DerivedParameter(variable_name = "dy",
                          type_ = "double",
		      defining_expression_in_cpp = """(this->yhigh - this->ylow) / this->my""")
    parameters.append(dy)

    dz = DerivedParameter(variable_name = "dz",
                          type_ = "double",
		      defining_expression_in_cpp = """(this->zhigh - this->zlow) / this->mz""")
    parameters.append(dz)
    
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list = \
    _parameters_accessors_checks()

def starter_pac():

    from finess.params import finess, weno, dim3_rotated_1d, append_pac_from_module
    
    pac = ([], [], [])
    
    append_pac_from_module(pac, finess)
    append_pac_from_module(pac, weno)   
    append_pac_from_module(pac, dim3_rotated_1d)
       
    return pac 

