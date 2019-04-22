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
include("test_ccacopf.jl")           ## test Γ sensitivity calcs

## test sensitivity on case9
const case = "case9"
opfdata = load_case(case, path, other=false);
const nbus = length(opfdata.buses)
const ngen = length(opfdata.generators)
const nload = sum(opfdata.buses.bustype .== 1)
sm = OPF.s_acopf_model(opfdata)
sm = OPF.acopf_solve(sm, opfdata)
sm_eval = setup(sm.m);               ## stochastic model evaluator
sm_zbar = deepcopy(sm_eval.last_x);  ## stochastc model equilibrium z̄
OPF.acopf_outputAll(sm, opfdata)
include("test_sensitivities.jl")     ## test Γ sensitivity calcs



Y = computeAdmittanceMatrix(opfdata)
m_idx = OPF.model_idx(opfdata)
z = deepcopy(sm_zbar);      nz = length(z)
x = deepcopy(z[m_idx[:x]]); nx = length(x)
y = deepcopy(z[m_idx[:y]]); ny = length(y)
Pnet, Qnet = OPF.PQnet(sm, opfdata)
rP, rQ = PFE(x, z, Y, Pnet, Qnet, m_idx)
F = fill(NaN, nx)
J = fill(NaN, nx, nx)
PFE!(F, x, z, Y, Pnet, Qnet, m_idx); F
PFE_J!(J, x, z, Y, Pnet, Qnet, m_idx, opfdata.BusIdx, opfdata.BusGenerators); J
PFE_wrap!(F, x) = PFE!(F, x, z, Y, Pnet, Qnet, m_idx)
PFE_J_wrap!(J, x) = PFE_J!(J, x, z, Y, Pnet, Qnet, m_idx, opfdata.BusIdx, opfdata.BusGenerators)
sol = nlsolve(PFE_wrap!, PFE_J_wrap!, x, iterations=5_000)
xx = deepcopy(sol.zero)

function Gamma_fd(x, z, Y, Pnet, Qnet, m_idx, epsilon)
    ## setup
    nz = length(z)
    nx = length(x)
    ny = length(m_idx[:y])
    nd = length(m_idx[:d])
    nbus = length(Pnet); @assert(nbus == length(Qnet)); @assert(nbus == size(Y,1))
    #### output
    Gamma_FD = fill(NaN, nx, nd)
    xs = fill(NaN, nx, nd)

    ## compute
    for i in 1:nd
        ## overwrite
        x_ = deepcopy(x)
        z_ = deepcopy(z)
        Pnet_ = deepcopy(Pnet)
        Qnet_ = deepcopy(Qnet)

        ## perturb
        if i <= nbus
            Pnet_[i] -= epsilon
        elseif i > nbus
            Qnet_[mod1(i, nbus)] -= epsilon
        end
        println(i)
        println("Pnet diff", Pnet-Pnet_)
        println("Qnet diff", Qnet-Qnet_)

        ## solve
        PFE_wrap!(F, x) = PFE!(F, x, z_, Y, Pnet_, Qnet_, m_idx)
        PFE_J_wrap!(J, x) = PFE_J!(J, x, z_, Y, Pnet_, Qnet_, m_idx, opfdata.BusIdx, opfdata.BusGenerators)
        sol = nlsolve(PFE_wrap!, PFE_J_wrap!, x_, iterations=5_000, ftol=epsilon/10)
        @info sol.f_converged == true
        xs[:,i] = sol.zero
        # println((xs[:,i] - x) / epsilon)
        println("NaNs: ", sum(isnan.((xs[:,i] - x) ./ epsilon)))
        println()
        Gamma_FD[:,i] = (xs[:,i] - x) / epsilon
    end
    return Gamma_FD, xs
end

G_FD, xs_ = Gamma_fd(xx, z, Y, Pnet, Qnet, m_idx, 1e-4)
_,_, dF = jac_z_alg(z, Y, opfdata.BusIdx, opfdata.BusGenerators, m_idx, false)
G = dF[:dF_dx] \ -dF[:dF_dd]
G-G_FD

## test
@test norm(Gamma - Gamma_FD) <= 10*tol
[norm(G[:,i] - G_[:,i]) for i = 1:ny]