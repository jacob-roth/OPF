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

function fpacopf_solve(case_name::String, constr_limit_scales::Array{T,1}, ratelimits::Array{T,1},
                           options0::Dict, fileout::String) where T <: AbstractFloat
    options = deepcopy(options0)
    casedata = load_case(case_name, path; other=true)
    opfdata = casedata.opf
    physdata = casedata.phys

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
            for key in [:solvetime, :objvalue, :iters_total, :iters_feasibility_phase, :iters_optimization_phase, :num_line_failure_constraints_added, :time_wall_rate_eval, :time_wall_main_nlp, :time_cpu_main_nlp, :time_cpu_main_nlp_in_ipopt, :time_cpu_main_nlp_func_eval, :gencost]
                println(key, " --> ", optimal_values[key])
            end
            optimal_values[:rates] = opfmodel_exitrates.other[:rates]
            optimal_values[:kktrates] = opfmodel_exitrates.other[:kktrates]
            optimal_values[:caputil] = opfmodel_exitrates.other[:caputil]
            optimal_values[:status] = string(opfmodel_exitrates.status)
            scale_fmt = @sprintf("%0.0d", Int(round(100(options[:constr_limit_scale]-1.0), digits=0)))
            rate_fmt  = @sprintf("%0.0e", options[:ratelimit])
            file_out = joinpath(fileout, "emergency=$(scale_fmt)pct", "case=$case_name", "$(rate_fmt)")
            mkpath(file_out)
            write_optimal_values(file_out * "/", optimal_values)
        end
    end
    return nothing
end

path = joinpath(dirname(@__FILE__), "..", "data")
fileout = joinpath(dirname(@__FILE__), "___temp___") #joinpath(dirname(@__FILE__), "output", "scalability_runs")

options = DefaultOptions()
options[:emergencylimit] = NaN
options[:ratelimit]      = NaN
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = true
options[:shed_load]      = false
options[:print_level]    = 5
options[:linear_solver]  = "ma27"
options[:allow_pf_infeas]= true
options[:VOLL]           = 1e+6
options[:tol]            = 1e-6

#
#------- NOTE -------#
# In pglib_opf_case300_ieee, line 179: from 1201 -> to 120 has negative reactance (-0.3697)
# N-0 ACOPF Objective value: 564008.4103087195. Generation cost = 564008.4071948407USD/hr
# if we delete the line: Objective value: 563640.9450685976. Generation cost = 563640.9419548993USD/hr
# if we make the reactance positive (+0.3697): Objective value: 563744.0633950386. Generation cost = 563744.0602812819USD/hr
# pglib_opf_case300_ieee_modified makes the reactance positive (+0.3697)
# all other casedata is unchanged
#
constr_limit_scales = [1.20]
acopf_maxrate = Dict(
    "pglib_opf_case118_ieee" => 7.713e-9, # actual: 7.757e-9. we use maximum(kktrates)
    "pglib_opf_case200_activ" => 6340.6, # actual: 4947.6. we use maximum(kktrates)
    "pglib_opf_case300_ieee_modified" => 6.754e-26, # actual: 6.044e-24. we use maximum(kktrates)
    "pglib_opf_case500_goc" => 374.0, # actual: 314.1. we use maximum(kktrates)
    "pglib_opf_case1354_pegase" => 310863.0, # actual: 221222.5, we use maximum(kktrates)
)
shortname = Dict(
    "pglib_opf_case118_ieee" => "118",
    "pglib_opf_case200_activ" => "200",
    "pglib_opf_case300_ieee_modified" => "300",
    "pglib_opf_case500_goc" => "500",
    "pglib_opf_case1354_pegase" => "1354",
)

RUN_CASES = ["pglib_opf_case118_ieee", "pglib_opf_case200_activ", "pglib_opf_case300_ieee_modified", "pglib_opf_case500_goc", "pglib_opf_case1354_pegase", "pglib_opf_case118_ieee"]

function run_scalability_runs()
    for case_name in RUN_CASES
        ratelimits = [Inf, 1e-1, 1e-2].*acopf_maxrate[case_name]
        fpacopf_solve(case_name, constr_limit_scales[1:1], ratelimits, options, fileout)
    end
end

function summarize_scalability_runs()
    for case_name in RUN_CASES
        ratelimits = [Inf, 1e-1, 1e-2].*acopf_maxrate[case_name]
        casedata = load_case(case_name, path; other=true)
        opfdata = casedata.opf
        physdata = casedata.phys
        total_Pd = sum(opfdata.buses.Pd)
        total_Qd = sum(opfdata.buses.Qd)
        gencost_acopf = Inf
        @printf("%6s", shortname[case_name]);
        for (idx, rate) in enumerate(ratelimits)
            scale_fmt = @sprintf("%0.0d", Int(round(100(options[:constr_limit_scale]-1.0), digits=0)))
            rate_fmt  = @sprintf("%0.0e", rate)
            outpath = joinpath(fileout, "emergency=$(scale_fmt)pct", "case=$case_name", "$(rate_fmt)")
            @assert readdlm("$outpath/status.csv")[1] == "Optimal"
            time_wall_main_nlp = readdlm("$outpath/time_wall_main_nlp.csv")[1]
            time_wall_rate_eval = readdlm("$outpath/time_wall_rate_eval.csv")[1]
            active_load_shed = sum(abs.(readdlm("$outpath/Ps.csv")))
            reactive_load_shed = sum(abs.(readdlm("$outpath/Qs.csv")))
            gencost = readdlm("$outpath/gencost.csv")[1]
            if isinf(rate)
                gencost_acopf = gencost
            end
            gencost_increase = 100.0*(gencost - gencost_acopf)/gencost_acopf
            load_shed = 100.0*(active_load_shed + reactive_load_shed)/(total_Pd + total_Qd)

            @printf(" & %8.1f", time_wall_main_nlp)
            if isinf(rate)
                continue
            end
            #@printf(" & %8.1f & %8.2f & %8.2f", time_wall_rate_eval, load_shed, gencost_increase)
            @printf(" & %8.1f & %8.2f", time_wall_rate_eval, load_shed)
        end
        @printf("\n")
    end
    nothing
end

function convert_pglib_case_to_mpccases(casepath::String)
    bus = readdlm(casepath * ".bus")
    gens = readdlm(casepath * ".gen")
    
    # the phys file
    nbus = size(bus, 1)
    ngen = size(gens, 1)
    phys = zeros(nbus, 4)
    phys[:,1] .= bus[:,1]
    phys[bus[:,2] .> 1,2] .= 0.0531
    phys[bus[:,2] .> 1,3] .= 0.05
    phys[bus[:,2] .== 1,3] .= 0.005
    phys[:,4] .= 0.01
    writedlm(casepath * ".phys", phys)

    # add additional columns to gen file
    num_new_cols = 21 - size(gens, 2)
    if num_new_cols > 0
        new_gens = zeros(ngen, 21)
        new_gens[:,1:size(gens,2)] .= gens
        writedlm(casepath * ".gen", new_gens)
    end
    nothing
end

run_scalability_runs()
summarize_scalability_runs()
