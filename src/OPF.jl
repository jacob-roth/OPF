# module OPF
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, NLsolve
using MatrixNetworks, StructArrays, Distributed ## for setting n-1 limits
using MPCCases
using Printf
using Pkg
## based on: https://github.com/StructJuMP/StructJuMP.jl/tree/master/examples/PowerGrid

mutable struct OPFModel
    m::JuMP.Model
    status::Symbol
    kind::Symbol
end
export OPFModel

include("default.jl")
export DefaultOptions, DefaultAdjustments

include("constraints.jl")
export add_p_constraint!, add_q_constraint!
export add_line_current_constraint!, add_line_power_constraint!

include("acpf_model.jl")
export acpf_model

include("acopf_model.jl")
export acopf_model

include("s_acopf_model.jl")
export s_acopf_model

include("acopf_zip_model.jl")
export acopf_zip_model

include("scacopf_model.jl")
export scacopf_model, get_operating_points

include("scacopf_n1_limits.jl")
export set_n1_limits

include("n1_limits.jl")
export get_n1_limits

include("util.jl")
export acopf_solve, acopf_initialPt_IPOPT
export acopf_outputAll, get_values
export RGL_id, RGL_idx, model_idx
export PQnet
export get_opfmodeldata
export get_contingencies

include("pfe.jl")
export P_i, Q_i
export PF, PFE, PFE!, PFE_J!
export PF_vec, PFE_vec, PFE_vec!

include("jacobian.jl")
export jac_z_num
export jac_z_alg_ew
export dStilde_dVtilde, jac_z_alg_vec
export jac_z

include("sensitivities.jl")
export get_Gamma, get_Gamma_fd, get_Gamma_ew

# end # module
