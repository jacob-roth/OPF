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
acopf_outputAll(opfmodel_acopf, opfdata, options) # Objective value: 125947.8781174833. Generation cost = 125947.87945987795USD/hr

##
## exitrates
##

## UNCONSTRAINED
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125949.81448816974. Generation cost = 125949.81587988648USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3 (line 58 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125969.23837945235. Generation cost = 125969.23974533858USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=25pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 5e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125976.71461805045. Generation cost = 125976.71598378317USD/hr
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
acopf_outputAll(opfmodel_acopf, opfdata, options) # Objective value: 125948.01275209706. Generation cost = 125948.01409450876USD/hr

##
## exitrates
##

## UNCONSTRAINED (line 125 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125948.01275120632. Generation cost = 125948.01411818918USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/ntap-nshunt/emergency=20pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1 (line 125 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125948.01275120632. Generation cost = 125948.01411818918USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/ntap-nshunt/emergency=20pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2 (line 113 & 125 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125952.14720358636. Generation cost = 125952.14859522428USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/ntap-nshunt/emergency=20pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125976.00455718602. Generation cost = 125976.00594848282USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/ntap-nshunt/emergency=20pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (line 3 & 37 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 5e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125984.65968474935. Generation cost = 125984.66107594887USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/ntap-nshunt/emergency=20pct/5e-4/118bus_5e-4"
write_optimal_values(file_out, optimal_values)

## 1e-5 (infeasible)
# options[:shed_load]      = true
# options[:ratelimit]      = 1e-5
# opfmodel = acopf_model(opfdata, options)
# opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
# acopf_outputAll(opfmodel_exitrates, opfdata, options) # 0
# optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
# file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/ntap-nshunt/emergency=20pct/1e-5/118bus_1e-5"
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
acopf_outputAll(opfmodel_acopf, opfdata, options) # Objective value: 125947.8781174833. Generation cost = 125947.87945987795USD/hr

##
## exitrates
##

## UNCONSTRAINED
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2 (line 68 negative) XXX-this case-XXX
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125953.47300979352. Generation cost = 125953.47440147224USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3 (line 68 negative) XXX-this case-XXX
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125978.57879647097. Generation cost = 125978.58016213935USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=15pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (line 182 infeasible; overall infeasible)
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
acopf_outputAll(opfmodel_acopf, opfdata, options) # Objective value: 125947.8781174833. Generation cost = 125947.87945987795USD/hr

##
## exitrates
##

## UNCONSTRAINED (line 37 negative, line 125 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1 (line 37 negative,  line 125 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2 (line 37 negative)
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125954.66211299953. Generation cost = 125954.66348103638USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3 (line 37 negative)
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125981.18010396104. Generation cost = 125981.18146957144USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (line 37 negative; infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 5e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: -1.8728181817365117e-6. Generation cost = 143250.68300539983USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=12-5pct/5e-4/118bus_5e-4"
write_optimal_values(file_out, optimal_values)

## -----------------------------------------------------------------------------
## get exitrate point @ 1.10x scaling
## -----------------------------------------------------------------------------

options = DefaultOptions()
options[:emergencylimit] = 1.10
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
acopf_outputAll(opfmodel_acopf, opfdata, options) # Objective value: 125947.8781174833. Generation cost = 125947.87945987795USD/hr

##
## exitrates
##

## UNCONSTRAINED (line 37 negative, line 125 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = Inf
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=10pct/inf/118bus_inf"
write_optimal_values(file_out, optimal_values)

## 1e-1 (line 37 negative,  line 125 infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 1e-1
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=10pct/1e-1/118bus_1e-1"
write_optimal_values(file_out, optimal_values)

## 1e-2
options[:shed_load]      = true
options[:ratelimit]      = 1e-2
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125955.87030524723. Generation cost = 125955.87167330092USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=10pct/1e-2/118bus_1e-2"
write_optimal_values(file_out, optimal_values)

## 1e-3 (line 68 negative)
options[:shed_load]      = true
options[:ratelimit]      = 1e-3
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125984.48608740517. Generation cost = 125984.48747893567USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=10pct/1e-3/118bus_1e-3"
write_optimal_values(file_out, optimal_values)

## 5e-4 (infeasible)
options[:shed_load]      = true
options[:ratelimit]      = 5e-4
opfmodel = acopf_model(opfdata, options)
opfmodel_exitrates = acopf_solve_exitrates(opfmodel, opfdata, options)
acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: -1.8728181817365117e-6. Generation cost = 143250.68300539983USD/hr
optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
file_out = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus-lowdamp/emergency=10pct/5e-4/118bus_5e-4"
write_optimal_values(file_out, optimal_values)
