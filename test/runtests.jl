## -----------------------------------------------------------------------------
## constants
## -----------------------------------------------------------------------------
const path = pwd() * "/cases/"
const tol = 1e-9
const plotting = false

## -----------------------------------------------------------------------------
## environment
## -----------------------------------------------------------------------------
import Pkg
Pkg.activate(dirname(dirname(dirname(path))))
Pkg.instantiate()
using Test
using MPCCases, Printf
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, Distributions
using NLsolve
include("../src/OPF.jl") # stand-in for "using OPF"

## -----------------------------------------------------------------------------
## test cases
## -----------------------------------------------------------------------------
for c in ["case9", "case30", "case118"]
    case = c
    opfdata = load_case(case, path, other=false);
    nbus = length(opfdata.buses);
    ngen = length(opfdata.generators);
    nload = sum(opfdata.buses.bustype .== 1);

    include("test_zip.jl")               ## test zip load model
    include("test_det-sto.jl")           ## compare deterministic and "stochastic"
    include("test_compare.jl")           ## compare with Anirudh's model
    include("test_pfe.jl")               ## test power flow equation calcs (vectorized, MatPower, and entrywise)
end

include("test_n1.jl")