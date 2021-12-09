using Distributed
@everywhere PROJ_DIR = joinpath(dirname(@__FILE__), "..")
@everywhere using Pkg
@everywhere Pkg.activate("$PROJ_DIR")
@everywhere Pkg.instantiate()
@everywhere begin
    using MPCCases, StructArrays, LinearAlgebra, ForwardDiff, Printf, SharedArrays, JuMP, Ipopt
    using TimerOutputs
    include("$PROJ_DIR/src/default.jl")
    include("$PROJ_DIR/src/exitrates.jl")
end
import Pkg; Pkg.activate("$PROJ_DIR"); Pkg.instantiate()
include("$PROJ_DIR/src/OPF.jl")
using DelimitedFiles

## -----------------------------------------------------------------------------
## data input
## -----------------------------------------------------------------------------

path = joinpath(dirname(@__FILE__), "..", "data")
fileout = joinpath(dirname(@__FILE__), "___temp___") #joinpath(dirname(@__FILE__), "output", "numerical_performance_118_tol=1e-06")
case_name = "mpc_lowdamp_pgliblimits"

options = DefaultOptions()
options[:emergencylimit] = NaN
options[:ratelimit]      = NaN
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = true
options[:shed_load]      = false
options[:print_level]    = 5
options[:allow_pf_infeas]= false
options[:linear_solver]  = "ma27"
options[:tol]            = 1e-6
#---------------------- NOTE ----------------------
# after review-1, the results in Table II were obtained
# using Ipoptv0.6.5 (see Project.toml [compat])
# with
#     em = Model(solver = IpoptSolver(print_level=options[:print_level]))
# inside `compute_exitrate_exact`
# and with
#     em = Model(solver = IpoptSolver(print_level=options[:print_level],max_iter=1000))
# inside `compute_exitrate_kkt`
# and with
# options[:tol] = 1e-9 for the main_nlp
# options[:linear_solver] = "mumps"
#--------------------------------------------------

function set_optimalvalues(case_name::String, constr_limit_scales::Array{T,1}, ratelimits::Array{T,1},
                           options0::Dict, fileout::String) where T <: AbstractFloat
    options = deepcopy(options0)
    casedata = load_case(case_name, path; other=true)
    opfdata = casedata.opf
    physdata = casedata.phys
    models = []

    ##
    ## exitrates
    ##

    for i in eachindex(constr_limit_scales)
        options[:constr_limit_scale] = constr_limit_scales[i]
        j_models = []
        for j in eachindex(ratelimits)
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
            file_out = joinpath(fileout, "emergency=$(scale_fmt)pct", "$(rate_fmt)")
            mkpath(file_out)
            write_optimal_values(file_out * "/", optimal_values)
        end
        push!(models, j_models)

        ## N-1
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
        file_out = joinpath(fileout, "N_1")
        mkpath(file_out)
        write_optimal_values(file_out * "/", optimal_values)
    end
    return models
end

function check_reliability(case_name::String, constr_limit_scales::Array, ratelimits::Array, options0::Dict)
    casedata = load_case(case_name, path; other=true)
    opfdata = casedata.opf
    options = deepcopy(options0)

    for i in eachindex(constr_limit_scales)
        options[:constr_limit_scale] = constr_limit_scales[i]
        for j in eachindex(ratelimits)
            options[:ratelimit] = ratelimits[j]
            scale_fmt = @sprintf("%0.0d", Int(round(100(options[:constr_limit_scale]-1.0), digits=0)))
            rate_fmt  = @sprintf("%0.0e", options[:ratelimit])
            file_out = fileout * "emergency=$(scale_fmt)pct/$(rate_fmt)"
            pg_ref = readdlm("$(file_out)/Pg.csv")[:,1]
            # pg_ref = readdlm("/home/asubramanyam/research/OPF/test/___temp___/N_1/Pg.csv")[:,1]
            feasible_ctgs = check_Nminus1_feasibility(pg_ref, opfdata, options)
            println("constr_limit_scale = $(options[:constr_limit_scale]), ratelimit = $(options[:ratelimit]) => reliability = ", length(feasible_ctgs), "/", length(opfdata.lines))
        end
    end
end

constr_limit_scales = [1.20]
ratelimits = [Inf, 1e-9, 1e-12, 1e-15]
models = set_optimalvalues(case_name, constr_limit_scales[1:1], ratelimits, options, fileout)
check_reliability(case_name, constr_limit_scales, ratelimits, options)
