"""Export parameter_list, accessor_list, check_list for two-fluid plasma equations.

These (free) parameters are the same identified as in ``Approximate Riemann
solver for the two-fluid plasma model'', U. Shumlak, J. Loverich, JCP 2003.
"""

from __future__ import absolute_import

def _parameters_accessors_checks():

    from finess.params import Parameter, DerivedParameter, Accessor,\
         Check, CheckGreaterEqual, CheckGreaterThan

    parameters = []
    checks     = []

    # Name of plasma model (for now, only g05 is implemented)
    model_name = Parameter(variable_name = "model_name",
                   section = "plasma",
    		       name    = "model_name",
    		       type_   = "std::string",
    		       default_value = "g05")
    parameters.append(model_name)		       

    # User supplied parameters that define electron-ion masses
#   ion_mass = Parameter(variable_name = "ion_mass",
#                     section = "plasma",
#                     name = "ion_mass",
#                     type_ = "double",
#                     default_value = 1.0
#                     )
#   parameters.append(ion_mass)
#   checks.append(CheckGreaterThan(ion_mass, 0.0))

    mass_ratio = Parameter(variable_name = "mass_ratio",
                      section = "plasma",
                      name = "mass_ratio",
                      type_ = "double",
                      default_value = 1.0
                      )
    parameters.append(mass_ratio)
    checks.append( CheckGreaterThan(mass_ratio, 1.0) )

    # Derived parameter: electron mass
#   electron_mass = DerivedParameter(variable_name = "elc_mass",
#                         type_    = "double",
#       defining_expression_in_cpp = """this->ion_mass/this->mass_ratio""")
#   parameters.append(electron_mass)


    debye_length = Parameter( variable_name = "debye_length",
                      section = "plasma",
                      name = "debye_length",
                      type_ = "double",
                      default_value = 1.0
                      )
    parameters.append(debye_length)
    checks.append( CheckGreaterThan(debye_length, 0.0) )

    larmor_radius = Parameter( variable_name = "larmor_radius",
                      section = "plasma",
                      name = "larmor_radius",
                      type_ = "double",
                      default_value = 1.0
                      )
    parameters.append(larmor_radius)
    checks.append( CheckGreaterThan(larmor_radius, 0.0) )


    # Parameters used for Perfectly Hyperbolic Maxwell (PHM)
    cc_factor = Parameter(variable_name = "cc_factor",
                      section = "plasma",
                      name = "cc_factor",
                      type_ = "double",
                      default_value = 1.1
                      )
    parameters.append( cc_factor )
    checks.append( CheckGreaterThan(cc_factor, 1.0) )

    # Derived parameter (should probably mix this with Maxwell, but eh ... )
    cc_sqd = DerivedParameter(variable_name = "cc_sqd",
                          type_    = "double",
        defining_expression_in_cpp = """this->cc_factor*this->cc_factor""")
    parameters.append(cc_sqd)


    clean_electric_field = Parameter(variable_name = "clean_E_field",
                   section = "plasma",
                   name    = "clean_E_field",
                   type_   = "int")
    parameters.append(clean_electric_field)
 
    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list =  \
    _parameters_accessors_checks()

