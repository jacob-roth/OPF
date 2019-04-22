## -----------------------------------------------------------------------------
## constants
## -----------------------------------------------------------------------------
const path = pwd() * "/cases/"
const case = "case9"
const tol = 1e-9

## -----------------------------------------------------------------------------
## environment
## -----------------------------------------------------------------------------
import Pkg
Pkg.activate(dirname(dirname(dirname(path))))
Pkg.instantiate()
using Test
using MPCCases, Printf, MAT
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, Distributions
using NLsolve
# using OPF
include("../src/OPF.jl")

## -----------------------------------------------------------------------------
## load
## -----------------------------------------------------------------------------
opfdata = load_case(case, path, other=false);
nbus = length(opfdata.buses)
ngen = length(opfdata.generators)
nload = sum(opfdata.buses.bustype .== 1)

## -----------------------------------------------------------------------------
## tests
## -----------------------------------------------------------------------------
include("test_det-sto.jl")           ## compare deterministic and "stochastic"
include("test_compare.jl")           ## compare with Anirudh's model
include("test_pfe.jl")               ## test power flow equation calcs (vectorized, MatPower, and entrywise)
include("test_jacobian.jl")          ## test jacobian calcs (algebraic v numerical)
include("test_sensitivities.jl")     ## test Γ sensitivity calcs
include("test_cc_acopf.jl")          ## test cc-acopf model