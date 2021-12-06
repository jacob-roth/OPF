path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files/"
tol = 1e-9
plotting = false
temp = 1e-4
damping = 1.0
case_name = "mpc_lowdamp_pgliblimits"
l = 50
include("/Users/jakeroth/git/OPF/src/OPF2PD.jl")

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

casedata = load_case(case_name, path; other=true)
opfdata  = casedata.opf
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
opfmodeldata = get_opfmodeldata(opfdata, options)
opfmodeldata[:Y] = imag.(opfmodeldata[:Y])
solution = get_optimal_values(opfmodel_acopf.m, opfmodeldata)

line_type = "PiModelLine"
phys = physDefault
opf2pd("/Users/jakeroth/git/PowerDynamics.jl/examples/118busnotapnoshunt/118notapnoshunt.json", solution, opfmodeldata, line_type, phys)
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118busnotapnoshunt/Pnet.csv", solution[:Pnet])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118busnotapnoshunt/Qnet.csv", solution[:Qnet])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118busnotapnoshunt/Vm.csv", solution[:Vm])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118busnotapnoshunt/Va.csv", solution[:Va])
writedlm("/Users/jakeroth/git/PowerDynamics.jl/examples/118busnotapnoshunt/Y.csv", opfmodeldata[:Y])
