@testset "Vm sensitivities" begin
## parameters
epsilon = 1e-5
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

## algebraic
J_alg = OPF.jac_z(sm_zbar, data, :ew)
J_alg = OPF.jac_z(sm_zbar, data, :ew, false)
dFdy, dFdx = OPF.dFdy_dFdx_RGL(J_alg, J_RGL_idx, full)
# dFdy, dFdx = dFdy[y_RGL_idx, y_RGL_idx], dFdx[y_RGL_idx, :]
dFdy, dFdx = dFdy[1:nload, 1:nload], dFdx[1:nload, :]
Gamma = Matrix(dFdy) \ -Matrix(dFdx)

## finite difference
F = fill(NaN, length(y_RGL))
# Vtilde, Itilde, Stilde = OPF.PF(y_RGL, data)
# OPF.PFE_RGL!(F, y_RGL, data); F
# PFE_RGL_wrap!(F, y) = OPF.PFE_RGL!(F, y, data)
# sol = nlsolve(PFE_RGL_wrap!, y_RGL, autodiff = :forward, method=:newton, ftol=1e-5, iterations=10_000)
# sol = nlsolve(PFE_RGL_wrap!, y_RGL, method=:newton, ftol=1e-6, iterations=10_000)

OPF.PFE_RGL_real!(F, y_RGL, data); F
PFE_RGL_real_wrap!(F, y) = OPF.PFE_RGL_real!(F, y, data)
sol = nlsolve(PFE_RGL_real_wrap!, y_RGL, ftol=1e-8)
OPF.PFE_RGL_real!(F, sol.zero, data); F

if norm(sol.zero - y_RGL) <= epsilon
    y_RGL = sol.zero  ## "refine" OPFModel solution...
end
Gamma_FD = fill(NaN, ny, nx)
y_RGLs = fill(NaN, ny, nx)
data_ = deepcopy(data)

## compute
for i = 1:nx
    ## perturb
    if i <= nbus
        ii = i
        data_[:Pnet][ii] = data[:Pnet][ii] - epsilon
    elseif i > nbus
        ii = mod1(i, nbus)
        data_[:Qnet][ii] = data[:Qnet][ii] - epsilon
    end
    println("Pnet diff: ", data_[:Pnet] - data[:Pnet])
    println("Qnet diff: ", data_[:Qnet] - data[:Qnet])
    println("")

    ## solve
    # PFE_RGL_wrap!(F, y) = OPF.PFE_RGL!(F, y, data_)
    # sol = nlsolve(PFE_RGL_wrap!, y_RGL, method=:newton, ftol=1e-6, iterations=10_000)
    PFE_RGL_real_wrap!(F, y) = OPF.PFE_RGL_real!(F, y, data_)
    sol = nlsolve(PFE_RGL_real_wrap!, y_RGL, method=:newton)
    @assert(sol.f_converged == true)
    y_RGLs[:,i] = sol.zero
    Gamma_FD[:,i] = (y_RGLs[:,i] - y_RGL) / epsilon

    ## reset
    if i <= nbus
        ii = i
        data_[:Pnet][ii] = data[:Pnet][ii]
    elseif i > nbus
        ii = mod1(i, nbus)
        data_[:Qnet][ii] = data[:Qnet][ii]
    end
end

## test
@test norm(Gamma - Gamma_FD) <= 10*tol
[norm(Gamma[:,i] - Gamma_FD[:,i]) for i = 1:nx]
end # testset





G =  [1 2 3 4 5;
      6 7 8 9 10]
S =  [0 0 0    0    0;
      0 0 0    0    0;
      0 0 0.9  0.3  0.7;
      0 0 0.1  1    2.1;
      0 0 0.2  30   4]
SS =     [0.9  0.3  0.7;
          0.1  1    2.1  ;
          0.2  30   4  ]
GG = [3 4 5;
      8 9 10]
G * S * G'
GG * SS * GG'