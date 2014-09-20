import finess.params.finess
import finess.params.weno
import finess.params.dim2
import finess.params.mhd

from finess.params import generate_header_cpp

parameters = []
accessors = []
checks = []


#section [finess]
#to add later
#mrestart
#nrestart
#datafmt
parameters += finess.params.finess.parameter_list
accessors += finess.params.finess.accessor_list
checks += finess.params.finess.check_list


#section [weno]
parameters += finess.params.weno.parameter_list
accessors += finess.params.weno.accessor_list
checks += finess.params.weno.check_list


#section [grid]

parameters += finess.params.dim2.parameter_list
accessors += finess.params.dim2.accessor_list
checks += finess.params.dim2.check_list


#section [mhd]
parameters += finess.params.mhd.parameter_list
accessors += finess.params.mhd.accessor_list
checks += finess.params.mhd.check_list

#parameters specific to current app
import this_app_params
parameters += this_app_params.parameter_list
accessors += this_app_params.accessor_list
checks += this_app_params.check_list


header_filename, header_code, cpp_filename, cpp_code = generate_header_cpp(parameters, accessors, checks)


#print header_code

#print cpp_code



with open(header_filename, 'w') as f:
    f.write(header_code)
with open(cpp_filename, 'w') as f:
    f.write(cpp_code)





