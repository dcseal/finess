"""Export parameter_list, accessor_list, check_list for Maxwell's equations.

Currently, only speed of light: maxwell.cs_light.

This value can be set in an [maxwell] section of your parameters.ini code.

Sorry Maxwell, we did not capitalize your name!  Please accept this
forgiveness.

Derived parameters:

    cs_light_squared - speed of light squared.

"""

from __future__ import absolute_import

def _parameters_accessors_checks():

    from finess.params import Parameter, DerivedParameter, Accessor,\
         Check, CheckGreaterEqual, CheckGreaterThan

    parameters = []
    checks = []
    cs_light = Parameter(variable_name = "cs_light",
                      section       = "maxwell",
                      name          = "cs_light",
                      type_         = "double",
                      default_value = 1.0
            )
    parameters.append( cs_light )
    checks.append(CheckGreaterThan(cs_light, 0.0))

    # Derived parameters (performed here to cut back on FLOPs)
    cs_light_squared = DerivedParameter(variable_name = "cs_light_squared",
                          type_         = "double",
        defining_expression_in_cpp = """this->cs_light*this->cs_light""")
    parameters.append(cs_light_squared)

    return parameters, map(Accessor, parameters), checks

parameter_list, accessor_list, check_list =  \
    _parameters_accessors_checks()


