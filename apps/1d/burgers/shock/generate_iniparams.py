from finess.params import append_pac_from_module, write_to_header_cpp
import finess.params.dim1

pac = finess.params.dim1.starter_pac()

parameter_list, accessor_list, check_list = pac

if __name__ == "__main__":
    write_to_header_cpp(*pac)

