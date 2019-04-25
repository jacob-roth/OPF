module OPF
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, NLsolve
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

include("acopf_model.jl")
export acopf_model

include("s_acopf_model.jl")
export s_acopf_model

include("util.jl")
export acopf_solve, acopf_initialPt_IPOPT
export acopf_outputAll, get_values
export RGL_id, RGL_idx, model_idx
export PQnet

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

end # module