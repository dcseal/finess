from finess.params import append_pac_from_module, write_to_header_cpp
import finess.params.dim3_rotated_1d
import finess.params.mhd

pac = finess.params.dim3_rotated_1d.starter_pac()

#section [mhd]
append_pac_from_module(pac, finess.params.mhd)

parameter_list, accessor_list, check_list = pac

if __name__ == "__main__":
    write_to_header_cpp(*pac)

