const path = pwd() * "/cases/"
const tol = 1e-9
const plotting = false

import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl")
using DelimitedFiles

## -----------------------------------------------------------------------------
## get exitrate point @ 1.25x scaling
## -----------------------------------------------------------------------------

options = DefaultOptions()
options[:emergencylimit] = 1.25
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
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

## UNCONSTRAINED
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125962.44155523923USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3 (line 58 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125949.80489121478USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 5e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125976.70509064505USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/5e-4/118bus_5e-4"
write_optimal_values(file_out, optimal_values)

## 1e-5 (infeasible)
# options[:shed_load]      = true
# options[:ratelimit]      = 1e-5
# opfmodel = acopf_model(opfdata, options)
# opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
# acopf_outputAll(opfmodel_exitrates, opfdata, options)
# optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
# file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/1e-5/118bus_1e-5"
# write_optimal_values(file_out, optimal_values)

## -----------------------------------------------------------------------------
## get exitrate point @ 1.20x scaling
## -----------------------------------------------------------------------------

options = DefaultOptions()
options[:emergencylimit] = 1.2
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
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

## UNCONSTRAINED
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=20pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=20pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125951.39815208458USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=20pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125973.68660292463USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=20pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4
options[:shed_load]      = true
options[:ratelimit]      = 5e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125981.87468039014USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=20pct/5e-4/118bus_5e-4"
write_optimal_values(file_out, optimal_values)

## 1e-5 (infeasible)
# options[:shed_load]      = true
# options[:ratelimit]      = 1e-5
# opfmodel = acopf_model(opfdata, options)
# opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
# acopf_outputAll(opfmodel_exitrates, opfdata, options) # 0
# optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
# file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=20pct/1e-5/118bus_1e-5"
# write_optimal_values(file_out, optimal_values)

## -----------------------------------------------------------------------------
## get exitrate point @ 1.15x scaling
## -----------------------------------------------------------------------------

options = DefaultOptions()
options[:emergencylimit] = 1.15
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = true
options[:ratelimit]      = 1.0
options[:print_level]    = 1
opfdata = load_case("case118-n1-lowdamp", path; other=false)

## acopf
options[:shed_load]      = false
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
acopf_outputAll(opfmodel_acopf, opfdata, options) ## cost = 125947.8781174833USD/hr

##
## exitrates
##

## UNCONSTRAINED (line 58 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1 (line 58 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2 (line 58 infeaasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125953.46325314212USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3 (line 58 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125978.56920721379USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (infeasible)
# options[:shed_load]      = true
# options[:ratelimit]      = 5e-4
# opfmodel = acopf_model(opfdata, options)
# opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
# acopf_outputAll(opfmodel_exitrates, opfdata, options) # 0
# optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
# file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/5e-4/118bus_5e-4"
# write_optimal_values(file_out, optimal_values)

## 1e-5 (infeasible)
# options[:shed_load]      = true
# options[:ratelimit]      = 1e-5
# opfmodel = acopf_model(opfdata, options)
# opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
# acopf_outputAll(opfmodel_exitrates, opfdata, options)
# optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
# file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/1e-5/118bus_1e-5"
# write_optimal_values(file_out, optimal_values)

## -----------------------------------------------------------------------------
## get exitrate point @ 1.125x scaling
## -----------------------------------------------------------------------------

options = DefaultOptions()
options[:emergencylimit] = 1.125
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = true
options[:ratelimit]      = 1.0
options[:print_level]    = 1
opfdata = load_case("case118-n1-lowdamp", path; other=false)

## acopf
options[:shed_load]      = false
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
acopf_outputAll(opfmodel_acopf, opfdata, options) ## cost = 125947.8781174833USD/hr

##
## exitrates
##

## UNCONSTRAINED (line 58 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125947.86870886046USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125954.65227896454USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 125981.17050187051USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (infeasible: sheds load)
options[:shed_load]      = true
options[:ratelimit]      = 5e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # 0
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/5e-4/118bus_5e-4"
write_optimal_values(file_out, optimal_values)

## 1e-5 (infeasible)
# options[:shed_load]      = true
# options[:ratelimit]      = 1e-5
# opfmodel = acopf_model(opfdata, options)
# opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
# acopf_outputAll(opfmodel_exitrates, opfdata, options)
# optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
# file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/1e-5/118bus_1e-5"
# write_optimal_values(file_out, optimal_values)
