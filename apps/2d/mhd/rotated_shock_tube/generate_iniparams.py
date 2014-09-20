import finess.params.finess
import finess.params.weno
import finess.params.dim2
import finess.params.mhd


parameter_list = []
accessor_list = []
check_list = []


#section [finess]
#to add later
#mrestart
#nrestart
#datafmt
parameter_list += finess.params.finess.parameter_list
accessor_list += finess.params.finess.accessor_list
check_list += finess.params.finess.check_list


#section [weno]
parameter_list += finess.params.weno.parameter_list
accessor_list += finess.params.weno.accessor_list
check_list += finess.params.weno.check_list


#section [grid]

parameter_list += finess.params.dim2.parameter_list
accessor_list += finess.params.dim2.accessor_list
check_list += finess.params.dim2.check_list


#section [mhd]
parameter_list += finess.params.mhd.parameter_list
accessor_list += finess.params.mhd.accessor_list
check_list += finess.params.mhd.check_list

#parameters specific to current app
import this_app_params
parameter_list += this_app_params.parameter_list
accessor_list += this_app_params.accessor_list
check_list += this_app_params.check_list

if __name__ == "__main__":
    from finess.params import generate_header_cpp
    header_filename, header_code, cpp_filename, cpp_code = generate_header_cpp(parameter_list, accessor_list, check_list)
    
    with open(header_filename, 'w') as f:
        f.write(header_code)
    with open(cpp_filename, 'w') as f:
        f.write(cpp_code)





