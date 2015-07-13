from finess.params import append_pac_from_module, write_to_header_cpp
import finess.params.dim1
import finess.params.mhd
import finess.params.maxwell
import finess.params.plasma
#import this_app_params

pac = finess.params.dim1.starter_pac()

# section [mhd]
append_pac_from_module(pac, finess.params.mhd     )     # defines gamma

# section [maxwell]
append_pac_from_module(pac, finess.params.maxwell )     # defines speed of light

# section [plasma]
append_pac_from_module(pac, finess.params.plasma  )     # defines ion-mass ratio, and PHM modifiers

# parameters specific to current app (TODO - include the left/right values)
# append_pac_from_module(pac, this_app_params)

parameter_list, accessor_list, check_list = pac

if __name__ == "__main__":
    write_to_header_cpp(*pac)




