const path = "/Users/jakeroth/Desktop/planning-large-deviation/data/cases/118-files/"
const path = "/home/jroth/Projects/planning-large-deviation/data/cases/118-files/"
const tol = 1e-9
const plotting = false
const temp = 1e-4
const damping = 1.0
const case_name = "mpc_base_pgliblimits"
const case_name = "mpc_lowdamp_pgliblimits"
const l = 50

import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl"); using DelimitedFiles

## -----------------------------------------------------------------------------
## data input
## -----------------------------------------------------------------------------

options = DefaultOptions()
options[:emergencylimit] = NaN
options[:ratelimit]      = NaN
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = false
options[:shed_load]      = false
options[:print_level]    = 1
options[:constr_limit_scale] = 1.04
options[:pw_angle_limits]= false
options[:slack0]         = true

## -----------------------------------------------------------------------------
## 0 - load case
## -----------------------------------------------------------------------------

mkpath("load-case/")
casedata = load_case(case_name, path; other=true)
opfdata = load_case(case_name, path; other=false)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
acopf_outputAll(opfmodel_acopf, opfdata, options)
opfmodeldata = get_opfmodeldata(opfdata, options)
solution = get_optimal_values(opfmodel_acopf.m, opfmodeldata)
writedlm("load-case/Pnet0.csv", solution[:Pnet])
writedlm("load-case/Qnet0.csv", solution[:Qnet])
writedlm("load-case/Vm0.csv", solution[:Vm])
writedlm("load-case/Va0.csv", solution[:Va])
writedlm("load-case/grad_H_xbar0.csv", solution[:grad_H])
writedlm("load-case/hess_H_xbar0.csv", solution[:hess_H])
writedlm("load-case/Y0.csv", opfmodeldata[:Y])

## -----------------------------------------------------------------------------
## extract some data
## -----------------------------------------------------------------------------

const buses        = opfmodeldata[:buses]
const line         = opfmodeldata[:lines][l]
const nonLoadBuses = opfmodeldata[:nonLoadBuses]
const bus_ref      = opfmodeldata[:bus_ref]
const Y            = opfmodeldata[:Y]
# flowmax      = options[:constr_limit_scale]*(line.rateA/(abs(1.0/(line.x*im))*opfmodeldata[:baseMVA]))^2
const flowmax      = options[:constr_limit_scale]*(line.rateA^2 / (opfmodeldata[:baseMVA]^2 * abs2(1.0/(line.x*im))))
const nbus         = length(buses)
const nrow         = 2nbus - length(nonLoadBuses) - 1
const VMbar0       = solution[:Vm]
const VAbar0       = solution[:Va]
const baseMVA      = opfmodeldata[:baseMVA]

## -----------------------------------------------------------------------------
## 1 - solve energy minimization problem
## -----------------------------------------------------------------------------
mkpath("xbar/")
VMbar = VMbar0
VAbar = VAbar0

em = Model(solver = IpoptSolver(print_level=options[:print_level]))

#
# Variables
#
@variable(em, Vm[1:nbus])
@variable(em, -pi <= Va[1:nbus] <= pi)
setlowerbound(Va[bus_ref], VAbar[bus_ref])
setupperbound(Va[bus_ref], VAbar[bus_ref])
setlowerbound(Vm[bus_ref], VMbar[bus_ref])
setupperbound(Vm[bus_ref], VMbar[bus_ref])
for b in 1:nbus
    if b in nonLoadBuses
        setlowerbound(Vm[b], VMbar[b])
        setupperbound(Vm[b], VMbar[b])
    end
end

#
# Set objective to minimize energy function
#
@NLobjective(em, Min,
    -0.5*sum(Y[m,n]*Vm[m]*Vm[n]*cos(Va[m] - Va[n]) for m in 1:nbus for n in 1:nbus)
    - sum(solution[:Pnet][m]*Va[m] + solution[:Qnet][m]*log(Vm[m]) for m in 1:nbus))

#
# initial value
#
setvalue(getindex(em, :Vm), VMbar)
setvalue(getindex(em, :Va), VAbar)

#
# solve model
#
status = solve(em)

#
# Construct solution in reduced space
#
VMbar1 = getvalue(getindex(em,:Vm))
VAbar1 = getvalue(getindex(em,:Va))
VMbar_r = copy(VMbar1)
VAbar_r = copy(VAbar1)
for i in opfmodeldata[:nonLoadBuses][end:-1:1]
    splice!(VMbar_r, i)
end
splice!(VAbar_r, opfmodeldata[:bus_ref])

writedlm("xbar/Vm.csv", VMbar1)
writedlm("xbar/Va.csv", VAbar1)
writedlm("xbar/grad_H_xbar.csv", ∇H([VMbar_r; VAbar_r]; solution=solution, opfmodeldata=opfmodeldata))
writedlm("xbar/hess_H_xbar.csv", ∇2H_direct([VMbar  ; VAbar  ]; solution=solution, opfmodeldata=opfmodeldata))

## -----------------------------------------------------------------------------
## 2 - solve constrained energy minimization problem
## -----------------------------------------------------------------------------
mkpath("xstar/")
em = Model(solver = IpoptSolver(print_level=options[:print_level]))

#
# Variables
#
@variable(em, Vm[1:nbus])
@variable(em, -pi <= Va[1:nbus] <= pi)
setlowerbound(Va[bus_ref], VAbar[bus_ref])
setupperbound(Va[bus_ref], VAbar[bus_ref])
setlowerbound(Vm[bus_ref], VMbar[bus_ref])
setupperbound(Vm[bus_ref], VMbar[bus_ref])
for b in 1:nbus
    if b in nonLoadBuses
        setlowerbound(Vm[b], VMbar[b])
        setupperbound(Vm[b], VMbar[b])
    end
end

#
# Set objective to minimize energy function
#
@NLobjective(em, Min,
    -0.5*sum(Y[m,n]*Vm[m]*Vm[n]*cos(Va[m] - Va[n]) for m in 1:nbus for n in 1:nbus)
    - sum(solution[:Pnet][m]*Va[m] + solution[:Qnet][m]*log(Vm[m]) for m in 1:nbus))

#
# Line failure constraint
#
i = opfmodeldata[:BusIdx][line.from]
j = opfmodeldata[:BusIdx][line.to]
# @NLconstraint(em, fail, baseMVA*( flowmax - abs(Vm[i]^2 + Vm[j]^2 - (2*Vm[i]*Vm[j]*cos(Va[i]-Va[j]))) ) == 0)
@NLconstraint(em, fail, baseMVA*( flowmax - (Vm[i]^2 + Vm[j]^2 - (2*Vm[i]*Vm[j]*cos(Va[i]-Va[j]))) ) == 0)

#
# initial value
#
setvalue(getindex(em, :Vm), VMbar)
setvalue(getindex(em, :Va), VAbar)

#
# solve model
#
status = solve(em)

#
# Construct solution in reduced space
#
VMstar = getvalue(getindex(em,:Vm))
VAstar = getvalue(getindex(em,:Va))

VMstar_r = copy(VMstar)
VAstar_r = copy(VAstar)
for i in opfmodeldata[:nonLoadBuses][end:-1:1]
    splice!(VMstar_r, i)
end
splice!(VAstar_r, opfmodeldata[:bus_ref])

#
# exit point
#
exit_point = Dict()
exit_point[:Vm]  = VMstar
exit_point[:Va]  = VAstar
exit_point[:Vmr] = VMstar_r
exit_point[:Var] = VAstar_r
exit_point[:Pnet]= solution[:Pnet]
exit_point[:Qnet]= solution[:Qnet]

H_xstar      =          H([VMstar_r; VAstar_r]; solution=exit_point, opfmodeldata=opfmodeldata)
grad_H_xstar =         ∇H([VMstar_r; VAstar_r]; solution=exit_point, opfmodeldata=opfmodeldata)
hess_H_xstar = ∇2H_direct([VMstar  ; VAstar  ]; solution=exit_point, opfmodeldata=opfmodeldata)
grad_h_xstar =         ∇h([VMstar_r; VAstar_r]; solution=exit_point, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)
hess_h_xstar =        ∇2h([VMstar_r; VAstar_r]; solution=exit_point, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)
hess_h_xstar = 0.5 * (hess_h_xstar + hess_h_xstar') # Make the Hessian symmetric

Vmistar = exit_point[:Vm][i]
Vmjstar = exit_point[:Vm][j]
Vaistar = exit_point[:Va][i]
Vajstar = exit_point[:Va][j]
g1 = ( Vmistar-(Vmjstar*cos(Vaistar - Vajstar)))^2
g2 = ( Vmjstar-(Vmistar*cos(Vaistar - Vajstar)))^2
g3 = ( Vmistar* Vmjstar*sin(Vaistar - Vajstar))^2
g4 = (-Vmistar* Vmjstar*sin(Vaistar - Vajstar))^2
Kstar = norm(grad_H_xstar) / norm(grad_h_xstar)
gSg = sum([(a > 0 : S[a]: 0);
           (b > 0 : S[b]: 0);
           (c > 0 : S[c]: 0);
           (d > 0 : S[d]: 0)] .* [g1; g2; g3; g4]) * Kstar^2

writedlm("xstar/Vm.csv", VMstar)
writedlm("xstar/Va.csv", VAstar)
writedlm("xstar/grad_H_xstar.csv", grad_H_xstar)
writedlm("xstar/hess_H_xstar.csv", hess_H_xstar)
writedlm("xstar/grad_theta_xstar.csv", grad_h_xstar)
writedlm("xstar/hess_theta_xstar.csv", hess_h_xstar)
writedlm("xstar/gSg.csv", gSg)

# # calculate rate
# Kstar      = norm(grad_H_xstar) / norm(grad_h_xstar)
# tempM      = hess_H_xstar - (Kstar * hess_h_xstar)
# z1         = det(xbar[:hess_H])/det(tempM)
# z2         = grad_H_xstar' * (tempM\grad_H_xstar)
# z3         = sqrt(abs(z1/z2))
# kappa      = Kstar^3/z3
# prefactor  = z3 * norm(grad_H_xstar)^2 * (options[:damping]/sqrt(2*pi*options[:temperature]) )
# energydiff = H_xstar - xbar[:H_xbar]
# expterm    = max(exp(-energydiff/options[:temperature]), eps(0.0))
# caputil    = 100.0 * options[:constr_limit_scale] * sqrt(VMbar[i]^2 + VMbar[j]^2 - 2*VMbar[i]*VMbar[j]*cos(VAbar[i]-VAbar[j]))/sqrt(flowmax)
# if isnan(kappa)
#     (options[:print_level] >= 1) && println("warning: unable to compute EXACT exit rate for line ", l, " = (", i, ",", j, "). kappa = ", kappa, " (NaN)")
#     return nothing
# end
#
# ##
# ##
# exit_point[:K]          = Kstar
# exit_point[:prefactor]  = prefactor
# exit_point[:energydiff] = energydiff
# exit_point[:expterm]    = expterm
# exit_point[:kappa]      = kappa
# exit_point[:caputil]    = caputil
#
#
#
# line = opfmodeldata[:lines][l]
# flowmax = (options[:constr_limit_scale]*line.rateA/(abs(1.0/(line.x*im))*opfmodeldata[:baseMVA]))^2
# flowmax_JR = options[:constr_limit_scale]*(line.rateA^2 / (opfmodeldata[:baseMVA]^2 * abs2(1.0/(line.x*im))))
#
# solution = get_optimal_values(opfmodel_acopf.m, opfmodeldata)
# ep2 = compute_exitrate_exact(l, solution, opfmodeldata, options)  ## NOTE: JR - USE THIS ONE
# writedlm("Va10.csv", ep2[:Va])
# writedlm("Vm10.csv", ep2[:Vm])
# writedlm("Pnet10.csv", ep2[:Pnet])
# writedlm("Qnet10.csv", ep2[:Qnet])
# writedlm("Y10.csv", opfmodeldata[:Y])
