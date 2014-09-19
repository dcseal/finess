
from finess.params import Parameter, DerivedParameter, Accessor, Check, \
                          CheckGreaterEqual, CheckGreaterThan, \
			  CheckOneOf, EnumParameterType, generate_header_cpp


class_name = "IniParams"
global_variable_name = "global_ini_params"

parameters = []
accessors = []
checks = []


#section [finess]
ndims = Parameter(variable_name = "ndims",
                  section = "finess",
                  name = "ndims",
                  type_ = "int")
parameters.append(ndims)
checks.append(CheckOneOf(ndims, [1, 2, 3]))

nout = Parameter(variable_name = "nout",
                 section = "finess",
                 name = "nout",
                 type_ = "int",
                 default_value = 1)
parameters.append(nout)
checks.append(CheckGreaterThan(nout, 0))

tfinal = Parameter(variable_name = "tfinal",
                   section = "finess",
                   name = "tfinal",
                   type_ = "double")
parameters.append(tfinal)
checks.append(CheckGreaterEqual(tfinal, 0.0))

initial_dt = Parameter(variable_name = "initial_dt",
                       section = "finess",
                       name = "initial_dt",
                       type_ = "double")
parameters.append(initial_dt)
checks.append(CheckGreaterThan(initial_dt, 0.0))

max_dt = Parameter(variable_name = "max_dt",
                   section = "finess",
                   name = "max_dt",
                   type_ = "double")
parameters.append(max_dt)
checks.append(Check(cpp_code = """if(!(this->max_dt >= this->initial_dt))
    terminate("finess.max_dt should > finess.initial_dt");
"""))

desired_cfl = Parameter(variable_name = "desired_cfl",
                        section = "finess",
                        name = "desired_cfl",
                        type_ = "double")
parameters.append(desired_cfl)
checks.append(CheckGreaterThan(desired_cfl, 0.0))

max_cfl = Parameter(variable_name = "max_cfl",
                    section = "finess",
                    name = "max_cfl",
                    type_ = "double")
parameters.append(max_cfl)
checks.append(Check(cpp_code = """if(!(this->max_cfl >= this->desired_cfl))
    terminate("finess.max_cfl should > this->desired_cfl");
"""))

#max number of time steps per call to FinSolvexxx
nv = Parameter(variable_name = "nv",
               section = "finess",
               name = "nv",
               type_ = "int")
parameters.append(nv)
checks.append(CheckGreaterThan(nv, 0))


#TIGHTLY COUPLED with other C++ code
#Need to change the other C++ code (RunFinpack)
time_stepping_method = Parameter(variable_name = "time_stepping_method",
                                 section = "finess",
                                 name = "time_stepping_method",
                                 type_ = EnumParameterType(enum_scope_name = "TimeSteppingMethod",
                                             string_enumerator_dict = \
                                                 {"Runge-Kutta": "RK",
                                                  "Lax-Wendroff": "LxW",
                                                  "Multiderivative": "MD",
                                                  "User-Defined": "USER_DEFINED"},
                                             class_name = class_name))
parameters.append(time_stepping_method)

space_order = Parameter(variable_name = "space_order",
                        section = "finess",
                        name = "space_order",
                        type_ = "int")
parameters.append(space_order)
checks.append(CheckOneOf(space_order, [1, 3, 5, 7, 9, 11]))

time_order = Parameter(variable_name = "time_order",
                       section = "finess",
                       name = "time_order",
                       type_ = "int")
parameters.append(time_order)
checks.append(CheckOneOf(time_order, [1, 2, 3, 4, 5]))

#use_limiter removed

verbosity = Parameter(variable_name = "verbosity",
                      section = "finess",
                      name = "verbosity",
                      type_ = "int")
parameters.append(verbosity)
checks.append(CheckOneOf(verbosity, [0, 1]))

mcapa = Parameter(variable_name = "mcapa",
                  section = "finess",
                  name = "mcapa",
                  type_ = "int")
parameters.append(mcapa)
checks.append(CheckGreaterEqual(mcapa, 0))

maux = Parameter(variable_name = "maux",
                 section = "finess",
                 name = "maux",
                 type_ = "int")
parameters.append(maux)
checks.append(Check(cpp_code="""if(!(this->maux >= this->mcapa))
    terminate("finess.maux should >= finess.mcapa");
"""))

source_term = Parameter(variable_name = "source_term",
                        section = "finess",
                        name = "source_term",
                        type_ = "bool")
parameters.append(source_term)

meqn = Parameter(variable_name = "meqn",
                 section = "finess",
                 name = "meqn",
                 type_ = "int")
parameters.append(meqn)
checks.append(CheckGreaterEqual(meqn, 1))

#to add later
#mrestart
#nrestart
#datafmt


#section [grid]
mx = Parameter(variable_name = "mx",
               section = "grid",
               name = "mx",
               type_ = "int")
parameters.append(mx)
checks.append(CheckGreaterThan(mx, 0))

my = Parameter(variable_name = "my",
               section = "grid",
               name = "my",
               type_ = "int")
parameters.append(my)
checks.append(CheckGreaterThan(my, 0))

mbc = Parameter(variable_name = "mbc",
                section = "grid",
                name = "mbc",
                type_ = "int")
parameters.append(mbc)
checks.append(CheckGreaterEqual(mbc, 0))
#May need some additional checks for the coherence of e.g. space order, with mbc

xlow = Parameter(variable_name = "xlow",
                 section = "grid",
		 name = "xlow",
		 type_ = "double")
parameters.append(xlow)

xhigh = Parameter(variable_name = "xhigh",
                  section = "grid",
		  name = "xhigh",
		  type_ = "double")
parameters.append(xhigh)
checks.append(Check(cpp_code = """if(!(xhigh > xlow))
    terminate("grid.xhigh should > grid.xlow.");"""))


ylow = Parameter(variable_name = "ylow",
                 section = "grid",
		 name = "ylow",
		 type_ = "double")
parameters.append(ylow)

yhigh = Parameter(variable_name = "yhigh",
                  section = "grid",
		  name = "yhigh",
		  type_ = "double")
parameters.append(yhigh)
checks.append(Check(cpp_code = """if(!(yhigh > ylow))
    terminate("grid.yhigh should > grid.ylow.");"""))

dx = DerivedParameter(variable_name = "dx",
                      type_ = "double",
		      defining_expression_in_cpp = """(this->xhigh - this->xlow) / this->mx""")
parameters.append(dx)

dy = DerivedParameter(variable_name = "dy",
                      type_ = "double",
		      defining_expression_in_cpp = """(this->yhigh - this->ylow) / this->my""")
parameters.append(dy)



#section [mhd]
gamma = Parameter(variable_name = "gamma",
                  section = "mhd",
                  name = "gamma",
                  type_ = "double")
parameters.append(gamma)
checks.append(CheckGreaterThan(gamma, 0.0))

#section [initial]
rhol = Parameter(variable_name = "rhol",
                 section = "initial",
                 name = "rhol",
                 type_ = "double")
parameters.append(rhol)
checks.append(CheckGreaterThan(rhol, 0.0))

unl = Parameter(variable_name = "unl",
                section = "initial",
                name = "unl",
                type_ = "double")
parameters.append(unl)

utl = Parameter(variable_name = "utl",
                section = "initial",
                name = "utl",
                type_ = "double")
parameters.append(utl)

u3l = Parameter(variable_name = "u3l",
                section = "initial",
                name = "u3l",
                type_ = "double")
parameters.append(u3l)

pl = Parameter(variable_name = "pl",
               section = "initial",
               name = "pl",
               type_ = "double")
parameters.append(pl)
checks.append(CheckGreaterThan(pl, 0.0))

Bnl = Parameter(variable_name = "Bnl",
                section = "initial",
                name = "Bnl",
                type_ = "double")
parameters.append(Bnl)

Btl = Parameter(variable_name = "Btl",
                section = "initial",
                name = "Btl",
                type_ = "double")
parameters.append(Btl)

B3l = Parameter(variable_name = "B3l",
                section = "initial",
                name = "B3l",
                type_ = "double")
parameters.append(B3l)

rhor = Parameter(variable_name = "rhor",
                 section = "initial",
                 name = "rhor",
                 type_ = "double")
parameters.append(rhor)
checks.append(CheckGreaterThan(rhor, 0.0))

unr = Parameter(variable_name = "unr",
                section = "initial",
                name = "unr",
                type_ = "double")
parameters.append(unr)

utr = Parameter(variable_name = "utr",
                section = "initial",
                name = "utr",
                type_ = "double")
parameters.append(utr)

u3r = Parameter(variable_name = "u3r",
                section = "initial",
                name = "u3r",
                type_ = "double")
parameters.append(u3r)

pr = Parameter(variable_name = "pr",
               section = "initial",
               name = "pr",
               type_ = "double")
parameters.append(pr)
checks.append(CheckGreaterThan(pr, 0.0))

Bnr = Parameter(variable_name = "Bnr",
                section = "initial",
                name = "Bnr",
                type_ = "double")
parameters.append(Bnr)

Btr = Parameter(variable_name = "Btr",
                section = "initial",
                name = "Btr",
                type_ = "double")
parameters.append(Btr)

B3r = Parameter(variable_name = "B3r",
                section = "initial",
                name = "B3r",
                type_ = "double")
parameters.append(B3r)

angle = Parameter(variable_name = "angle",
                  section = "initial",
                  name = "angle",
                  type_ = "double")
parameters.append(angle)

#section [weno]

weno_version = Parameter(variable_name = "weno_version",
                   section = "weno",
		   name = "weno_version",
		   type_ = EnumParameterType( \
		             enum_scope_name = "WenoVersion",
			     string_enumerator_dict = \
			         {"FD": "FD",
				  "JS": "JS",
				  "Z": "Z"},
		             class_name = class_name),
	           default_value = "JS")
parameters.append(weno_version)

power_param = Parameter(variable_name = "power_param",
                       section = "weno",
		       name = "power_param",
		       type_ = "double",
		       default_value = 2.0)
parameters.append(power_param)



alpha_scaling = Parameter(variable_name = "alpha_scaling",
                          section = "weno",
                          name = "alpha_scaling",
                          type_ = "double",
                          default_value = 1.1)
parameters.append(alpha_scaling)
checks.append(CheckGreaterEqual(alpha_scaling, 1.0))


epsilon = Parameter(variable_name = "epsilon",
                    section = "weno",
                    name = "epsilon",
                    type_ = "double",
                    default_value = 1e-6)
parameters.append(epsilon)
checks.append(CheckGreaterThan(epsilon, 0.0))



#Something from section [finess]
global_alpha = Parameter(variable_name = "global_alpha",
                         section = "finess",
                         name = "global_alpha",
                         type_ = "bool",
                         default_value = False)
parameters.append(global_alpha)

output_dir = Parameter(variable_name = "output_dir",
                       section = "finess",
		       name = "output_dir",
		       type_ = "std::string",
		       default_value = "output")
parameters.append(output_dir)		       

accessors = map(Accessor, parameters)

header_filename, header_code, cpp_filename, cpp_code =     generate_header_cpp(class_name, global_variable_name, parameters, accessors, checks)


#print header_code

#print cpp_code



with open(header_filename, 'w') as f:
    f.write(header_code)
with open(cpp_filename, 'w') as f:
    f.write(cpp_code)





