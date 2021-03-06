const path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/30-files/"
const tol = 1e-9
const plotting = false

import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl")
using DelimitedFiles

## -----------------------------------------------------------------------------
## data input
## -----------------------------------------------------------------------------

constr_limit_scales = [1.10, 1.125, 1.15, 1.20, 1.25].^2
ratelimits = [Inf; exp10.(collect(-1:-1:-3)); 5e-4]
fileout = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/30bus-lowdamp/"
case_name = "mpc-lowdamp"

options = DefaultOptions()
options[:emergencylimit] = NaN
options[:ratelimit]      = NaN
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = true
options[:print_level]    = 1

function set_optimalvalues(case_name::String, constr_limit_scales::Array{T,1}, ratelimits::Array{T,1},
                           options0::Dict, fileout::String) where T <: AbstractFloat
    options = deepcopy(options0)
    casedata = load_case(case_name, path; other=true)
    opfdata = casedata.opf
    physdata = casedata.phys
    models = []

    ## acopf
    options[:shed_load] = false
    opfmodel = acopf_model(opfdata, options)
    opfmodel_acopf = acopf_solve(opfmodel, opfdata)
    acopf_outputAll(opfmodel_acopf, opfdata, options)
    # 118-n1-lowdamp: Objective value: 125947.8781174833. Generation cost = 125947.87945987795USD/hr
    # 30-mpc-lowdamp: Objective value: 604.249424491365. Generation cost = 604.249424491365USD/hr

    ##
    ## exitrates
    ##

    for i in eachindex(constr_limit_scales)
        options[:constr_limit_scale] = constr_limit_scales[i]
        j_models = []
        for j in eachindex(ratelimits)
            options[:shed_load] = true
            options[:ratelimit] = ratelimits[j]
            println("\nline scale = $(options[:constr_limit_scale])x")
            println("rate limit = $(@sprintf("%0.0e", options[:ratelimit]))\n")

            opfmodel = acopf_model(opfdata, options)
            opfmodel_exitrates = acopf_solve_exitrates(opfmodel, casedata, options)
            push!(j_models, opfmodel_exitrates)

            # rates = zeros(3, 41)
            # m = opfmodel_exitrates
            # for i in 1:41
            #     rates[:,i] .= [m.other[:prefactors][i]; m.other[:expterms][i]; m.other[:rates][i]]
            # end

            acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
            optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
            optimal_values[:rates] = opfmodel_exitrates.other[:rates]
            optimal_values[:status] = string(opfmodel_exitrates.status)
            scale_fmt = @sprintf("%0.0d", Int(round(100(options[:constr_limit_scale]-1.0), digits=0)))
            rate_fmt  = @sprintf("%0.0e", options[:ratelimit])
            file_out = fileout * "emergency=$(scale_fmt)pct/$(rate_fmt)/"
            mkpath(file_out)
            write_optimal_values(file_out, optimal_values)
        end
        push!(models, j_models)
    end
    return models
end

constr_limit_scales = round.(collect(range(1.20, stop=1.4, step=0.05)).^2, digits=2)
ratelimits = [Inf; exp10.(collect(range(3, stop=-4, step=-1)))]
models = set_optimalvalues(case_name, constr_limit_scales[1:1], ratelimits, options, fileout)

# set_optimalvalues(case_name, constr_limit_scales[1:end], ratelimits[1:end], options, fileout)
sumrates = zeros(length(models[1]))
for i in eachindex(models[1])
    sumrates[i] = sum(models[1][i].other[:rates])
end
