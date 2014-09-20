from __future__ import absolute_import

def read_params(parameter_filename, parameters_list):
    from ConfigParser import ConfigParser
    from finess.params import DerivedParameter
    config = ConfigParser()
    config.read(parameter_filename)
    d = dict(sum([[((s, opt), config.get(s, opt)) for opt in config.options(s)] for s in config.sections()], []))

    for p in parameters_list:
        if isinstance(p, DerivedParameter) or not (p.section, p.name) in d:
            continue
        if p.type_.type_string == 'int':
            d[p.section, p.name] = int(d[p.section, p.name])
        if p.type_.type_string == 'bool':
            d[p.section, p.name] = d[p.section, p.name] == 'true'
        if p.type_.type_string == 'double':
            d[p.section, p.name] = float(d[p.section, p.name])

    return d

    


