module OPF
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, Distributions
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

include("cc_acopf_model.jl")
export cc_acopf_model

include("s_acopf_model.jl")
export s_acopf_model

include("util.jl")
export acopf_solve, acopf_initialPt_IPOPT
export cc_acopf_solve
export acopf_outputAll, get_values
export RGL_id, RGL_idx, model_idx
export om_z_idx
export om_x_RGL_idx, om_y_RGL_idx, om_f_RGL_idx, om_pfe_RGL_idx, om_jac_RGL_idx
export PQnet

include("jacobian.jl")
export jac_z_num, jac_z_alg_ew, dStilde_dVtilde, jac_z_alg_vec
export jac_z_alg
export jac_z
export dFdy_dFdx_RGL

include("pfe.jl")
export PF, PFE_RGL!
export PF_real, PFE_RGL_real!

end # module