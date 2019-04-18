using JuMP
using Ipopt
using MathProgBase
using ForwardDiff
using LinearAlgebra
using Test

## 1/2 xAx + bx
## -----------------------------------------------------------------------------
A = [10.0 1.0; 1.0 9.0]
b = [2.0; 3.0]
tol = 1e-4
m = Model(solver=IpoptSolver())
@variable(m, x[1:2])
@variable(m, zeta[1:2])
setvalue(x, [0.0; 0.0])
@NLobjective(m, Min, sum(zeta[i]^2 for i = 1:2))
@NLconstraint(m, constraint, 0.5 * sum(A[i, j] * x[i] * x[j] for i in eachindex(x)
                                                       for j in eachindex(x))
                                     + sum(b[i] * x[i] for i in eachindex(x)) <= 10)
v = exp.(im .* x)
@NLconstraint(m, constraint_imag, sum(real.(abs(v[i])) for i = 1:2) <= 100)
## AD constraint
g = fill(NaN, 2)
d = JuMP.NLPEvaluator(m)
MathProgBase.initialize(d, [:Jac])
J = spzeros(1,2)
ij_jac = MathProgBase.jac_structure(d)
j = ones(length(ij_jac[1]))
MathProgBase.eval_jac_g(d, j, getvalue(x))
for elem in zip(ij_jac[1], ij_jac[1], j)
    J[elem[1], elem[2]] = elem[3]
end
for i = 1:2
    @NLconstraint(m, J[1,i] == zeta[i])
end

## solve and test
solve(m)
@test norm(getvalue(m[:x]) + A\b) <= tol
