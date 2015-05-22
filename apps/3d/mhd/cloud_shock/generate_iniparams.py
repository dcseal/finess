from finess.params import append_pac_from_module, write_to_header_cpp
import finess.params.dim3
import finess.params.mhd
import finess.params.artificial_viscosity
import this_app_params

pac = finess.params.dim3.starter_pac()

#section [mhd]
append_pac_from_module(pac, finess.params.mhd)

#section [av]
append_pac_from_module(pac, finess.params.artificial_viscosity)

append_pac_from_module(pac, this_app_params)

parameter_list, accessor_list, check_list = pac

if __name__ == "__main__":
    write_to_header_cpp(*pac)

