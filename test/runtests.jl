## -----------------------------------------------------------------------------
## constants
## -----------------------------------------------------------------------------
const path = pwd() * "/cases/"
const case = "case30"
# const case = "case5"
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
using OPF
# include("../src/OPF.jl")

## -----------------------------------------------------------------------------
## load
## -----------------------------------------------------------------------------
opfdata = load_case(case, path, other=false);

## -----------------------------------------------------------------------------
## compare deterministic and "stochastic"
## -----------------------------------------------------------------------------
## deterministic acopf
dm = OPF.acopf_model(opfdata)
dm = OPF.acopf_solve(dm, opfdata)
dm_eval = setup(dm.m);
dxbar = deepcopy(dm_eval.last_x);
OPF.acopf_outputAll(dm, opfdata)

## stochastic acopf
sm = OPF.sacopf_model(opfdata)
sm = OPF.acopf_solve(sm, opfdata)
sm_eval = setup(sm.m);
sxbar = deepcopy(sm_eval.last_x);
OPF.acopf_outputAll(sm, opfdata)

@testset "indexing" begin
@test norm([ getvalue(dm.m[:Pg]);
             getvalue(dm.m[:Qg]);
             getvalue(dm.m[:Vm]);
             getvalue(dm.m[:Va]) ] - dxbar) <= tol
@test norm([ getvalue(sm.m[:Pg]);
             getvalue(sm.m[:Qg]);
             getvalue(sm.m[:Vm]);
             getvalue(sm.m[:Va]);
             getvalue(sm.m[:Pd]);
             getvalue(sm.m[:Qd]) ] - sxbar) <= tol
@test norm([ getvalue(getindex(dm.m, :Pg));
             getvalue(getindex(dm.m, :Qg));
             getvalue(getindex(dm.m, :Vm));
             getvalue(getindex(dm.m, :Va)) ] - dxbar) <= tol
@test norm([ getvalue(getindex(sm.m, :Pg));
             getvalue(getindex(sm.m, :Qg));
             getvalue(getindex(sm.m, :Vm));
             getvalue(getindex(sm.m, :Va));
             getvalue(getindex(sm.m, :Pd));
             getvalue(getindex(sm.m, :Qd)) ] - sxbar) <= tol
end

## -----------------------------------------------------------------------------
## compare with other model
## -----------------------------------------------------------------------------
include("compare.jl")

## -----------------------------------------------------------------------------
## get Jacobian
## -----------------------------------------------------------------------------
## setup
Y = computeAdmittanceMatrix(opfdata)
idx_full = get_idx_sets(opfdata, full=true)
idx_red = get_idx_sets(opfdata, full=false)
@testset "get Jacobian" begin
for idx in [idx_full, idx_red]
## numerical
J_numerical = jac_x(sxbar, model=sm_eval)
J_numerical = J_numerical[idx[:f], :]

## algebraic (compare to MatPower)
Vtilde, Itilde, Stilde = PF(sxbar, Y, idx)
SVm, SVa = dStilde_dV(Vtilde, Y)
SVa_ = matread("cases/SVa.mat")
SVm_ = matread("cases/SVm.mat")
@test norm(SVm - SVm_["SVm"]) <= tol
@test norm(SVa - SVa_["SVa"]) <= tol
J_algebraic = jac_x(sxbar, Y, idx)
J_algebraic = J_algebraic[idx[:f], :]

## test
@test norm(J_algebraic - J_numerical) <= tol
end
end # testset

## -----------------------------------------------------------------------------
## get sensitivities
## -----------------------------------------------------------------------------
@testset "consistency between algebraic/numerical"
for idx in [idx_full, idx_red]
dFdy_n, dFdx_n = dFdy_dFdx(sxbar, idx, model=sm_eval)
dydx_numerical = Matrix(dFdy_n) \ Matrix(dFdx_n)
display(cond(Matrix(dFdy_n)))

dFdy_a, dFdx_a = dFdy_dFdx(sxbar, Y, idx)
dydx_algebraic = Matrix(dFdy_a) \ Matrix(dFdx_a)
display(cond(Matrix(dFdy_a)))
@test norm(dydx_algebraic - dydx_numerical) <= tol
end
end # testset (alg-num)

@testset "consistency between full/reduced" begin
## full
dFdy_f, dFdx_f = dFdy_dFdx(sxbar, Y, idx_full)
dydx_f = Matrix(dFdy_f) \ Matrix(dFdx_f)
display(cond(Matrix(dFdy_f)))

## reduced
dFdy_r, dFdx_r = dFdy_dFdx(sxbar, Y, idx_red)
dydx_r = Matrix(dFdy_r) \ Matrix(dFdx_r)
display(cond(Matrix(dFdy_r)))

rem = filter(x -> x ∉ idx_red[:f], idx_full[:f])
dFdy_f[7:end, ]

end # testset (full-reduced)




## -----------------------------------------------------------------------------
## test model
## -----------------------------------------------------------------------------
YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmittances(opfdata)
admittance = Dict()
admittance[:YffR] = YffR; admittance[:YffI] = YffI
admittance[:YttR] = YttR; admittance[:YttI] = YttI
admittance[:YftR] = YftR; admittance[:YftI] = YftI
admittance[:YtfR] = YtfR; admittance[:YtfI] = YtfI
admittance[:YshR] = YshR; admittance[:YshI] = YshI






## -----------------------------------------------------------------------------
## test sensitivities
## -----------------------------------------------------------------------------
Y = computeAdmittanceMatrix(opfdata)
G = real.(Y)
B = imag.(Y)
bus_idx, gen_idx = RGL_idx(opfdata)
bus_RGL = [bus_idx[:R]; bus_idx[:G]; bus_idx[:L]]
bus_RGL_inv = sortperm(bus_RGL)
G_RGL = G[bus_RGL, bus_RGL]
B_RGL = B[bus_RGL, bus_RGL]
vars, pars, eqns = jac_idx(opfdata)

## test PF equations consistency
Vtilde, Itilde, Stilde = PF_real(sxbar, G_RGL, B_RGL, pars)
Vtilde_, Itilde_, Stilde_ = PF(sxbar, Y, idx)
@test norm(real.(Vtilde_[bus_RGL]) - Vtilde.R) <= tol; @test norm(real.(Vtilde_) - Vtilde.R[bus_RGL_inv]) <= tol
@test norm(imag.(Vtilde_[bus_RGL]) - Vtilde.I) <= tol; @test norm(imag.(Vtilde_) - Vtilde.I[bus_RGL_inv]) <= tol
@test norm(real.(Itilde_[bus_RGL]) - Itilde.R) <= tol; @test norm(real.(Itilde_) - Itilde.R[bus_RGL_inv]) <= tol
@test norm(imag.(Itilde_[bus_RGL]) - Itilde.I) <= tol; @test norm(imag.(Itilde_) - Itilde.I[bus_RGL_inv]) <= tol
@test norm(real.(Stilde_[bus_RGL]) - Stilde.R) <= tol; @test norm(real.(Stilde_) - Stilde.R[bus_RGL_inv]) <= tol
@test norm(imag.(Stilde_[bus_RGL]) - Stilde.I) <= tol; @test norm(imag.(Stilde_) - Stilde.I[bus_RGL_inv]) <= tol

## test PF jacobian consistency
dStilde_dVm, dStilde_dVa = dStilde_dVtilde_real(sxbar, G_RGL, B_RGL, pars)
dStilde_dVm_, dStilde_dVa_ = dStilde_dVtilde(Vtilde_, Y)
@test norm(real.(dStilde_dVm_[bus_RGL, bus_RGL]) - dStilde_dVm.R) <= tol
@test norm(imag.(dStilde_dVm_[bus_RGL, bus_RGL]) - dStilde_dVm.I) <= tol
@test norm(real.(dStilde_dVa_[bus_RGL, bus_RGL]) - dStilde_dVa.R) <= tol
@test norm(imag.(dStilde_dVa_[bus_RGL, bus_RGL]) - dStilde_dVa.I) <= tol



const eps = 1e-7
idx = get_idx_sets(opfdata, full=false)
Pd = opfdata.buses.Pd
Qd = opfdata.buses.Qd
dydx_fd = zeros(length(idx[:y]), length(idx[:x]))
# for i in eachindex(Pd)

Pnet, Qnet = PQnet(sm, opfdata)
data = Dict()
data[:FF] = fill(NaN, 2length(opfdata.buses))
data[:idx] = idx_full
data[:Y] = deepcopy(Y)
data[:xmodel] = deepcopy(sxbar)
data[:Pnet] = deepcopy(Pnet)
data[:Qnet] = deepcopy(Qnet)
y = deepcopy(sxbar[idx[:y]])

F = fill(NaN, length(y))
PFE!(F, y, data)
PFE_wrap!(F, y) = PFE!(F, y, data)
sol = nlsolve(PFE_wrap!, y)