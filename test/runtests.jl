## -----------------------------------------------------------------------------
## constants
## -----------------------------------------------------------------------------
const path = pwd() * "/cases/"
const case = "case30"
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
using SparseArrays, LinearAlgebra
using NLsolve
# using OPF
include("../src/OPF.jl")

## -----------------------------------------------------------------------------
## load
## -----------------------------------------------------------------------------
opfdata = load_case(case, path, other=false);
const nbus = length(opfdata.buses)
const ngen = length(opfdata.generators)

## -----------------------------------------------------------------------------
## tests
## -----------------------------------------------------------------------------
include("test_det-sto.jl")           ## compare deterministic and "stochastic"
include("test_compare.jl")           ## compare with Anirudh's model
include("test_jacobian.jl")          ## test jacobian calcs (algebraic v numerical)
include("test_pfe.jl")               ## test power flow equation calcs (vectorized, MatPower, and entrywise)
# include("test_sensitivities.jl")     ## test Γ sensitivity calcs

## test sensitivity on case9
const case = "case9"
opfdata = load_case(case, path, other=false);
const nbus = length(opfdata.buses)
const ngen = length(opfdata.generators)
sm = OPF.sacopf_model(opfdata)
sm = OPF.acopf_solve(sm, opfdata)
sm_eval = setup(sm.m);               ## stochastic model evaluator
sm_zbar = deepcopy(sm_eval.last_x);  ## stochastc model equilibrium z̄
OPF.acopf_outputAll(sm, opfdata)
include("test_sensitivities.jl")     ## test Γ sensitivity calcs
