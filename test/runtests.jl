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
const nload = sum(opfdata.buses.bustype .== 1)

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
const nload = sum(opfdata.buses.bustype .== 1)
sm = OPF.sacopf_model(opfdata)
sm = OPF.acopf_solve(sm, opfdata)
sm_eval = setup(sm.m);               ## stochastic model evaluator
sm_zbar = deepcopy(sm_eval.last_x);  ## stochastc model equilibrium z̄
OPF.acopf_outputAll(sm, opfdata)
include("test_sensitivities.jl")     ## test Γ sensitivity calcs

## -----------------------------------------------------------------------------
## ccopf
## -----------------------------------------------------------------------------
full = false

## setup
Y = computeAdmittanceMatrix(opfdata)
z_idx = OPF.om_z_idx(opfdata)
J_RGL_idx = OPF.om_jac_RGL_idx(opfdata)
x_RGL_idx = OPF.om_x_RGL_idx(J_RGL_idx);                nx = length(x_RGL_idx)
y_RGL_idx = OPF.om_y_RGL_idx(J_RGL_idx, full)[1:nload]; ny = length(y_RGL_idx)
f_RGL_idx = OPF.om_f_RGL_idx(J_RGL_idx, full)[1:nload]; nf = length(f_RGL_idx)
y_RGL = sm_zbar[y_RGL_idx]
Pnet, Qnet = OPF.PQnet(sm, opfdata)

## data
data = Dict()
data[:Y]          = Y
data[:z_idx]      = z_idx
data[:y_RGL_idx]  = y_RGL_idx
data[:f_RGL_idx]  = f_RGL_idx
data[:opfmodel_z] = sm_zbar
data[:opfdata]    = opfdata
data[:Pnet]       = Pnet
data[:Qnet]       = Qnet
data[:Sigma]      = Matrix(Diagonal(ones(nx)))

## options
options = Dict()
