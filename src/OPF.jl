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
export acopf_solve
export acopf_outputAll
export acopf_initialPt_IPOPT
export get_idx_and_id
export filter_namedtuple
export nonunique

include("jacobian.jl")
export jac_x
export PF, dStilde_dV
export get_idx_sets
export dFdy_dFdx

end # module