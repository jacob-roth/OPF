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
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = true
options[:ratelimit]      = 1.0
options[:print_level]    = 1
opfdata = load_case("case30", path; other=false)

## acopf
# options[:shed_load]      = false
# opfmodel = acopf_model(opfdata, options)
# opfmodel_acopf = acopf_solve(opfmodel, opfdata)
# acopf_outputAll(opfmodel_acopf, opfdata, options)

##
## exitrates
##

## as close to "noshed" as possible
options[:shed_load]      = false
options[:ratelimit]      = 0.27325 ## this is close to the limit on line 29
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/30bus/noshed/30bus_noshed"
write_optimal_values(file_out, optimal_values)

## 1e-1
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/30bus/1e-1/30bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-3
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/30bus/1e-3/30bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 1e-5
options[:shed_load]      = true
options[:ratelimit]      = 1e-5
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options)
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/30bus/1e-5/30bus_1e-5"
write_optimal_values(file_out, optimal_values)