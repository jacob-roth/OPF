const path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files/"
# const path = "/home/jroth/Projects/planning-large-deviation/data/cases/118-files/"
const tol = 1e-9
const plotting = false
const temp = 1e-4
const damping = 1.0
# const case_name = "mpc_base_pgliblimits"
const case_name = "mpc_lowdamp_pgliblimits"
const l = 50

import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl"); using DelimitedFiles

## -----------------------------------------------------------------------------
## data input
## -----------------------------------------------------------------------------

options = DefaultOptions()
options[:emergencylimit] = NaN
options[:ratelimit]      = NaN
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = false
options[:print_level]    = 1
options[:constr_limit_scale] = 1.04
options[:pw_angle_limits]= false
options[:slack0]         = true

## -----------------------------------------------------------------------------
## 0 - load case
## -----------------------------------------------------------------------------

mkpath("load-case/")
casedata = load_case(case_name, path; other=true)
# opfdata = load_case(case_name, path; other=true)
opfdata  = casedata.opf
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
acopf_outputAll(opfmodel_acopf, opfdata, options)
opfmodeldata = get_opfmodeldata(opfdata, options)
opfmodeldata[:Y] = imag.(opfmodeldata[:Y])
solution = get_optimal_values(opfmodel_acopf.m, opfmodeldata)

line_type = "PiModelLine"
opf2pd("/Users/jakeroth/git/PowerDynamics.jl/examples/118bus/118notap.json", solution, opfmodeldata, line_type, phys)

writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118bus/Pnet.csv", solution[:Pnet])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118bus/Qnet.csv", solution[:Qnet])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118bus/Vm.csv", solution[:Vm])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118bus/Va.csv", solution[:Va])
# writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/pglib_118bus/grad_H_xbar.csv", solution[:grad_H])
# writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/pglib_118bus/hess_H_xbar.csv", solution[:hess_H])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118bus/Y.csv", opfmodeldata[:Y])
