import finess.params.dim2
import finess.params.mhd
import this_app_params

pac = finess.params.dim2.starter_pac()

#section [mhd]
finess.params.append_pac_from_module(pac, finess.params.mhd)

#parameters specific to current app
finess.params.append_pac_from_module(pac, this_app_params)

parameter_list, accessor_list, check_list = pac

if __name__ == "__main__":
    from finess.params import generate_header_cpp
    header_filename, header_code, cpp_filename, cpp_code = generate_header_cpp(*pac)
    
    with open(header_filename, 'w') as f:
        f.write(header_code)
    with open(cpp_filename, 'w') as f:
        f.write(cpp_code)





