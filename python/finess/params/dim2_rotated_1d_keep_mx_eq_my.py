""" Rotated 1D.  Keeping mx == my.
   
    Free parameters:
        Name in C++     Name in ini file
        mx              grid.mx
        angle           initial.angle   the angle the '1D direction' makes with x-axis
        mbc             grid.mbc
    
    Derived parameters:
        xlow    0
        xhigh   1 / cos(angle)
        ylow    0
        yhigh   1 / sin(angle)
        my      mx
        dx      (xhigh - xlow) / mx
        dy      (yhigh - ylow) / my
    
    Reference:
        doi:10.1006/jcph.2000.6519      Section 6.3.1
"""


from __future__ import absolute_import

def _parameters_accessors_checks():
    from finess.params import Parameter, DerivedParameter, Accessor,\
         Check, CheckGreaterEqual, CheckGreaterThan, CheckLessThan
    from math import pi

    parameters = []
    checks = []
    mx = Parameter(variable_name = "mx",
                   section = "grid",
                   name = "mx",
                   type_ = "int")
    parameters.append(mx)
    checks.append(CheckGreaterThan(mx, 0))

    angle = Parameter(variable_name = 'angle',
                      section = 'initial',
                      name = 'angle',
                      type_ = 'double')
    parameters.append(angle)
    checks.append(CheckGreaterThan(angle, 0))
    checks.append(CheckLessThan(angle, pi / 2))


    mbc = Parameter(variable_name = "mbc",
                    section = "grid",
                    name = "mbc",
                    type_ = "int")
    parameters.append(mbc)
    checks.append(CheckGreaterEqual(mbc, 0))

   
    my = DerivedParameter(variable_name = "my",
                          type_ = "double",
		                  defining_expression_in_cpp = """this->mx""")
    parameters.append(my)
    checks.append(CheckGreaterThan(my, 0))
    
    xlow = DerivedParameter(variable_name = "xlow",
                            type_ = "double",
                            defining_expression_in_cpp = "0.0")
    parameters.append(xlow)

    xhigh = DerivedParameter(variable_name = "xhigh",
                             type_ = "double",
                             defining_expression_in_cpp = "1 / std::cos(this->angle)")
    parameters.append(xhigh)

    ylow = DerivedParameter(variable_name = "ylow",
                            type_ = "double",
                            defining_expression_in_cpp = "0.0")
    parameters.append(ylow)

    yhigh = DerivedParameter(variable_name = "yhigh",
                             type_ = "double",
                             defining_expression_in_cpp = "1 / std::sin(this->angle)")
    parameters.append(yhigh)
    
   
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


def starter_pac():
    from finess.params import finess, weno, dim2_rotated_1d_keep_mx_eq_my, append_pac_from_module
    
    pac = ([], [], [])
    
    append_pac_from_module(pac, finess)
    append_pac_from_module(pac, weno)   
    append_pac_from_module(pac, dim2_rotated_1d_keep_mx_eq_my)
       
    return pac 

