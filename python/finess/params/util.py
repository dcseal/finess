from __future__ import absolute_import

def read_params(parameter_filename, parameters_list):
    from ConfigParser import ConfigParser
    from finess.params import DerivedParameter
    config = ConfigParser()
    config.read(parameter_filename)
    params = dict(sum([[((s, opt), config.get(s, opt)) for opt in config.options(s)] for s in config.sections()], []))

    for p in parameters_list:
        if isinstance(p, DerivedParameter) or not (p.section, p.name) in params:
            continue
        if p.type_.type_string == 'int':
            params[p.section, p.name] = int(params[p.section, p.name])
        if p.type_.type_string == 'bool':
            params[p.section, p.name] = params[p.section, p.name] == 'true'
        if p.type_.type_string == 'double':
            params[p.section, p.name] = float(params[p.section, p.name])

    if ("finess", "output_dir") in params and \
        params["finess", "output_dir"] != "":
        output_dir = params["finess", "output_dir"]
    else:
        output_dir = "output"
    params["finess", "output_dir"] = output_dir

   
    return params

    


