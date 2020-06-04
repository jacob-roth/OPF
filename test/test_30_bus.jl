const path = pwd() * "/cases/"
const tol = 1e-9
const plotting = false

import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl")
using DelimitedFiles


options = DefaultOptions()
options[:constr_limit_scale] = 1.44
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = true
options[:ratelimit]      = 1.0
options[:print_level]    = 1
options[:temperature]    = 1e-4
options[:damping]        = 1.0
casedata = load_case("case30", path; other=true)
opfdata  = casedata.opf
fileoutpath = ""

models = []
for ratelimit in [Inf, 1e3, 1e2, 1e1, 1e-0, 1e-1, 1e-2, 1e-3, 1e-4]
    options[:ratelimit] = ratelimit
        println("\n--------------------")
        println("line scale = $(options[:constr_limit_scale])x")
        println("rate limit = $(@sprintf("%0.0e", options[:ratelimit]))")
        println("--------------------\n")
    opfmodel = acopf_model(opfdata, options)
    opfmodel_exitrates = acopf_solve_exitrates(opfmodel, casedata, options)
    push!(models, opfmodel_exitrates)

    optimal_values = get_optimal_values(opfmodel_exitrates.m, get_opfmodeldata(opfdata, options))
    optimal_values[:rates] = opfmodel_exitrates.other[:rates]
    optimal_values[:status] = string(opfmodel_exitrates.status)
    scale_fmt = @sprintf("%0.0d", Int(round(100(options[:constr_limit_scale]-1.0), digits=0)))
    rate_fmt  = @sprintf("%0.0e", options[:ratelimit])
    file_out = fileoutpath * "emergency=$(scale_fmt)pct/$(rate_fmt)/"
    mkpath(file_out)
    write_optimal_values(file_out, optimal_values)
end

sumrates = zeros(length(models))
for i in eachindex(models)
    sumrates[i] = sum(models[i].other[:rates])
end

x1 = readdlm("/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/30bus-lowdamp/emergency=44pct/1e+02/rates.csv")
x2 = models[3].other[:rates]