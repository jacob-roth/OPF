const path = pwd() * "/cases/"
const tol = 1e-9
const plotting = false

import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl")

##
## get exitrate point
##

options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = false
options[:shed_load]      = true
options[:ratelimit]      = 1.0
options[:print_level]    = 1
opfdata = load_case("case118-n1-lowdamp", path; other=false)

## acopf
options[:shed_load]      = false
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
acopf_outputAll(opfmodel_acopf, opfdata, options) ## cost = 125_947.8781174833USD/hr

##
## exitrates
##

## as close to "noshed" as possible
options[:shed_load]      = false
# options[:ratelimit]      = 0.01 ## 4 lines
# options[:ratelimit]      = 0.0175 ## 4 lines
# options[:ratelimit]      = 0.01825 ## 4 lines
# options[:ratelimit]      = 0.0183 ## 4 lines
# options[:ratelimit]      = 0.0184 ## 4 lines
# options[:ratelimit]      = 0.018425 ## 4 lines
# options[:ratelimit]      = 0.01843 ## 4 lines
# options[:ratelimit]      = 0.018435 ## 4 lines
options[:ratelimit]      = 0.0184375 ## 4 lines
# options[:ratelimit]      = 0.018438 ## og
# options[:ratelimit]      = 0.018439 ## og
# options[:ratelimit]      = 0.01844 ## og
# options[:ratelimit]      = 0.01845 ## og
# options[:ratelimit]      = 0.0185 ## og
# options[:ratelimit]      = 0.019 ## og
# options[:ratelimit]      = 0.02 ## og
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/noshed/118bus_noshed"
write_optimal_values(file_out, optimal_values)

## 1e-1
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 1e-4 (infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/1e-4/118bus_1e-4"
write_optimal_values(file_out, optimal_values)

## 1e-5 (infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-5
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/1e-5/118bus_1e-5"
write_optimal_values(file_out, optimal_values)