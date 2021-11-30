module OPF
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, NLsolve, ForwardDiff
using MatrixNetworks, StructArrays, Distributed, SharedArrays ## for setting n-1 limits
using MPCCases
using Printf, DelimitedFiles
using Pkg
using Statistics
using TimerOutputs
## based on: https://github.com/StructJuMP/StructJuMP.jl/tree/master/examples/PowerGrid

mutable struct OPFModel
    m::JuMP.Model
    status::Symbol
    kind::Symbol
    other::Dict
end
export OPFModel

include("default.jl")
export DefaultOptions, DefaultAdjustments

include("constraints.jl")
export add_p_constraint!, add_q_constraint!
export add_line_current_constraint!, add_line_power_constraint!

include("acpf_model.jl")
export acpf_model

include("acopf_model_Pg.jl")
export acopf_model_Pg

include("acopf_model.jl")
export acopf_model

include("agg_branch_files.jl")
export write_new_limits, write_agg_branch_file

include("s_acopf_model.jl")
export s_acopf_model

include("acopf_zip_model.jl")
export acopf_zip_model

include("scacopf_model.jl")
export scacopf_model

include("scacopf_n1_limits.jl")
export set_n1_limits!, get_acopf_point, get_scacopf_point, get_contingency_point, update_limits!

include("util.jl")
export acopf_solve, acopf_solve_Pg, acopf_initialPt_IPOPT, acpf_solve
export acopf_outputAll
export scacopf_solve
export RGL_id, RGL_idx, model_idx
export PQnet
export get_opfmodeldata
export get_all_contingencies, subset_contingencies
export remove_line!, reinstate_line!
export get_flowmag2s, get_ratings
export get_islanding_buses, get_islanding_buses, get_nonislanding_contingencies
export get_point

include("pfe.jl")
export P_i, Q_i
export PF, PFE, PFE!, PFE_J!
export PF_vec, PFE_vec

include("jacobian.jl")
export jac_z_num
export jac_z_alg_ew
export dStilde_dVtilde, jac_z_alg_vec
export jac_z

include("sensitivities.jl")
export get_Gamma, get_Gamma_fd, get_Gamma_ew

include("exitrates.jl")
export acopf_solve_exitrates
export compute_exitrate_kkt, compute_exitrate_exact
export write_optimal_values, get_optimal_values
export compute_add_exitrates_parallel, compute_add_exitrates_serial
export compute_exitrate_exact_all_parallel, compute_exitrate_exact_all_serial

include("operatingdata.jl")
export get_operatingdata

include("OPF2PD.jl")
export opf2pd

# include("acopf_n1_limits.jl")
# export set_n1_limits!
# export get_nonislanding_lines
# export adjust_feas_ratings!, check_feasibility
# export adjust_solv_ratings!, check_solvability
# export remove_line!
# export reinstate_line!
# export get_opf_point
# export get_flowmag2s
# export get_ratings

end # module