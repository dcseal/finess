from finess.params import append_pac_from_module, write_to_header_cpp
import finess.params.dim2_rotated_1d_keep_mx_eq_my
import finess.params.mhd
import this_app_params

pac = finess.params.dim2_rotated_1d_keep_mx_eq_my.starter_pac()

#section [mhd]
append_pac_from_module(pac, finess.params.mhd)


parameter_list, accessor_list, check_list = pac

if __name__ == "__main__":
    write_to_header_cpp(*pac)




