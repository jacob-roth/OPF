module OPF
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays
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

include("opfmodel.jl")
export acopf_model

include("sopfmodel.jl")
export sacopf_model

include("util.jl")
export acopf_solve, acopf_initialPt_IPOPT
export acopf_outputAll, get_values
export RGL_id, RGL_idx
export om_z_idx
export om_x_RGL_idx, om_y_RGL_idx, om_pfe_RGL_idx, om_jac_RGL_idx
export PQnet

include("jacobian.jl")
export jac_x
export dFdy_dFdx_RGL
# export PF, PFE!, dStilde_dVtilde
end # module