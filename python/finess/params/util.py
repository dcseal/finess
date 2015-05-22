"""Provides functions related to params module, but used only in
plotting and not IniParams generation.
"""

from __future__ import absolute_import

def read_params(parameter_filename, parameters_list):
    """Read from parameter_filename.
    Interpret the contents as a .ini file.
    With the information from parameters_list, do type conversion if
    rhs of an assignment line is int/bool/double.
    Return a dict of type (str, str) --> int/bool/double/str. 
    parameters_list is the list that is used to generate IniParams.
    """

    from ConfigParser import ConfigParser
    from finess.params import DerivedParameter
    config = ConfigParser()
    config.read(parameter_filename)
    params = dict(sum([[((s, opt), config.get(s, opt)) for opt in config.options(s)] for s in config.sections()], []))

    for p in parameters_list:
        if isinstance(p, DerivedParameter) and p.dump_to_section == None:
            continue
        if not isinstance(p, DerivedParameter) and not (p.section, p.name) in params:
            continue
        if p.type_.type_string == 'int':
            params[p.section, p.name] = int(params[p.section, p.name])
        if p.type_.type_string == 'bool':
            params[p.section, p.name] = params[p.section, p.name] == 'true'
        if p.type_.type_string == 'double':
            params[p.section, p.name] = float(params[p.section, p.name])

    return params

    


