const path = "/home/jroth/Projects/planning-large-deviation/data/cases/118-files"
const path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files"
const path = "/home/asubramanyam/research/planning-large-deviation/data/cases/118-files"
const tol = 1e-9
const plotting = false

using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("..")
@everywhere Pkg.instantiate()
@everywhere begin
    using MPCCases, StructArrays, LinearAlgebra, ForwardDiff, Printf, SharedArrays, JuMP, Ipopt
    using TimerOutputs
    include("../src/default.jl")
    include("../src/exitrates.jl")
end
import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl")
using DelimitedFiles

## -----------------------------------------------------------------------------
## data input
## -----------------------------------------------------------------------------

fileout = "/home/jroth/Projects/planning-large-deviation/data/optimalvalues/118bus_lowdamp/"
fileout = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus_lowdamp/ytap_yshunt/emergency=4pct/"
fileout = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus_lowdamp/ytap_yshunt/"
fileout = "./___temp___/"
case_name = "mpc_lowdamp_pgliblimits"

options = DefaultOptions()
options[:emergencylimit] = NaN
options[:ratelimit]      = NaN
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = false
options[:shed_load]      = true
options[:print_level]    = 5

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

            acopf_outputAll(opfmodel_exitrates, opfdata, options) # Objective value: 125947.87811710747. Generation cost = 125947.87948380895USD/hr
            optimal_values = get_optimal_values(opfmodel_exitrates.m, opfmodel_exitrates.other[:opfmodeldata])
            optimal_values[:rates] = opfmodel_exitrates.other[:rates]
            optimal_values[:status] = string(opfmodel_exitrates.status)
            scale_fmt = @sprintf("%0.0d", Int(round(100(options[:constr_limit_scale]-1.0), digits=0)))
            rate_fmt  = @sprintf("%0.0e", options[:ratelimit])
            file_out = fileout * "emergency=$(scale_fmt)pct/$(rate_fmt)/"
            mkpath(file_out)
            write_optimal_values(file_out, optimal_values)
        end
        push!(models, j_models)

        ## N-1
        options[:ramp_pct] = 0.01
        casedata = load_case(case_name, path; other=true);
        opfmodeldata = get_opfmodeldata(casedata.opf, options)
        contingencies = get_all_contingencies(casedata.opf, options)
        scopfmodel = scacopf_model(casedata.opf, options, DefaultAdjustments(), contingencies)
        scopfmodel_acopf = scacopf_solve(scopfmodel, opfdata, options, contingencies)

        rates = zeros(length(opfmodeldata[:lines]))
        prefactors = zeros(length(opfmodeldata[:lines]))
        expterms = zeros(length(opfmodeldata[:lines]))
        opfmd = get_opfmodeldata(casedata, options)
        opfmd[:Y] = imag.(opfmd[:Y])
        optimal_values = get_optimal_values(scopfmodel_acopf.m, opfmd)
        pl = options[:print_level]
        options[:print_level] = 0
        for l in eachindex(opfmd[:lines])
            ep = compute_exitrate_exact(l, optimal_values, opfmd, options)
            if ep != nothing
                exitrate = ep[:prefactor] * ep[:expterm]
                rates[l] = exitrate
                prefactors[l] = ep[:prefactor]
                expterms[l] = ep[:expterm]
            end
        end
        options[:print_level] = pl
        optimal_values[:rates] = rates
        optimal_values[:status] = string(scopfmodel_acopf.status)
        file_out = fileout * "N_1/"
        mkpath(file_out)
        write_optimal_values(file_out, optimal_values)
    end
    return models
end

constr_limit_scales = [1.04] #round.(collect(range(1.0, stop=1.1, step=0.02)).^2, digits=2)
ratelimits = [Inf, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11] #[Inf; exp10.(collect(range(3, stop=-4, step=-1)))]
models1 = set_optimalvalues(case_name, constr_limit_scales[1:1], ratelimits, options, fileout)
ratelimits = [1e-12, 1e-16, 1e-20, 1e-24, 1e-28, 1e-30]
models2 = set_optimalvalues(case_name, constr_limit_scales[1:1], ratelimits, options, fileout)

ratelimits = [1e-30, 1e-40, 1e-50, 1e-60]
models = set_optimalvalues(case_name, constr_limit_scales[1:1], ratelimits, options, fileout)

# set_optimalvalues(case_name, constr_limit_scales[1:end], ratelimits[1:end], options, fileout)
sumrates = zeros(length(models[1]))
for i in eachindex(models[1])
    sumrates[i] = sum(models[1][i].other[:rates])
end


# file_out1 = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus_lowdamp/ytap_yshunt/emergency=4pct/1e-30/"
# file_out2 = "/Users/jakeroth/Desktop/planning-large-deviation/data/optimalvalues/118bus_lowdamp/ytap_yshunt/emergency=4pct/1e-30_v2/"
# x1 = readdlm(file_out1 * "rates.csv")
# x2 = readdlm(file_out2 * "rates.csv")
# norm(x1-x2)