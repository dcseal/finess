try:
    from finess.params import append_pac_from_module, write_to_header_cpp
except(ImportError):
    err  = '\nUnable to import module finess.params\n'
    err += 'Make sure that $FINESS/python/finess is part of the environment variable $PYTHONPATH\n'
    raise( ImportError(err) )

import finess.params.dim1

pac = finess.params.dim1.starter_pac()
parameter_list, accessor_list, check_list = pac

if __name__ == "__main__":
    write_to_header_cpp(*pac)

