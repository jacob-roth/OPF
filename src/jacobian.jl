"""
## `dFdy_dFdx_RGL`: partition jacobian `J` into `dFdy` and `dFdy` with `RGL`-ordering
### arguments:
    - `J::AbstractArray`: full jacobian (2nbus × (2ngen + 2nbus + 2nbus))
    - `jac_RGL_idx::Dict`: RGL index sets containing `f`, `x`, and `y`
### returns:
    - `dFdy::AbstractArray`: dF/dy
    - `dFdx::AbstractArray`: dF/dx
"""
function dFdy_dFdx_RGL(J::AbstractArray, jac_RGL_idx::Dict, full=false)
    if full
        fidx = [jac_RGL_idx[:f][:P_G]; jac_RGL_idx[:f][:P_L]; jac_RGL_idx[:f][:Q_RGL]]
        yidx = [jac_RGL_idx[:y][:Qg_RG]; jac_RGL_idx[:y][:Vm_L]; jac_RGL_idx[:y][:Va_GL]]
    else
        fidx = [jac_RGL_idx[:f][:P_G]; jac_RGL_idx[:f][:P_L]; jac_RGL_idx[:f][:Q_L]]
        yidx = [jac_RGL_idx[:y][:Vm_L]; jac_RGL_idx[:y][:Va_GL]]
    end
    xidx = [jac_RGL_idx[:x][:Pd_RGL]; jac_RGL_idx[:x][:Qd_RGL]]
    dFdy = J[fidx, yidx]
    dFdx = J[fidx, xidx]
    return dFdy, dFdx
end

## -----------------------------------------------------------------------------
## numerical
## -----------------------------------------------------------------------------
function jac_x(opfmodel_z::Array{Float64,1}; model::JuMP.NLPEvaluator)
    nconstr = MathProgBase.numconstr(model.m)
    nvar = MathProgBase.numvar(model.m)
    J = spzeros(nconstr, nvar)
    j!(J, opfmodel_z, model=model)
    return J
end

## -----------------------------------------------------------------------------
## algebraic (real, entrywise): http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/makeJac.html
## -----------------------------------------------------------------------------
function jac_x(opfmodel_z::AbstractArray, data::Dict=Dict(), full=true)
    """ based on: http://schevalier.com/wp-content/uploads/2017/02/Power-Flow-and-Covariance-Matrix.pdf """
    ## extract
    Y = data[:Y]
    G = real.(data[:Y])
    B = imag.(data[:Y])
    busIdx = data[:opfdata].BusIdx
    generators = data[:opfdata].generators
    z_idx = data[:z_idx]
    nbus = length(data[:opfdata].buses)
    ngen = length(data[:opfdata].generators)

    ## setup
    J = zeros(2nbus, 2nbus)
    Vm = opfmodel_z[z_idx[:Vm]]
    Va = opfmodel_z[z_idx[:Va]]

    ## compute
    for q = 1:nbus # P, Q; equations
        for b = 1:nbus # Vm, Va; buses
            h = busIdx[mod1(b, nbus)]
            k = busIdx[mod1(q, nbus)]
            if h == k
                IDX = [busIdx[x] for x in Y[h,:].nzind]
                P = Vm[h] * sum(Vm[kk] * ( G[h,kk] * cos(Va[h]-Va[kk]) + B[h,kk] * sin(Va[h]-Va[kk]) ) for kk in IDX)
                Q = Vm[h] * sum(Vm[kk] * ( G[h,kk] * sin(Va[h]-Va[kk]) - B[h,kk] * cos(Va[h]-Va[kk]) ) for kk in IDX)
                dP_dVa = -Q       - B[h,h] * Vm[h]^2
                dP_dVm =  P/Vm[h] + G[h,h] * Vm[h]
                dQ_dVa =  P       - G[h,h] * Vm[h]^2
                dQ_dVm =  Q/Vm[h] - B[h,h] * Vm[h]
            else
                dP_dVa =  Vm[h] * Vm[k] * ( G[h,k] * sin(Va[h]-Va[k]) - B[h,k] * cos(Va[h]-Va[k]) )
                dP_dVm =  Vm[h]         * ( G[h,k] * cos(Va[h]-Va[k]) + B[h,k] * sin(Va[h]-Va[k]) )
                dQ_dVa = -Vm[h] * Vm[k] * ( G[h,k] * cos(Va[h]-Va[k]) + B[h,k] * sin(Va[h]-Va[k]) )
                dQ_dVm =  Vm[h]         * ( G[h,k] * sin(Va[h]-Va[k]) - B[h,k] * cos(Va[h]-Va[k]) )
            end
            J[b, q]           = dP_dVm
            J[b, nbus+q]      = dP_dVa
            J[nbus+b, q]      = dQ_dVm
            J[nbus+b, nbus+q] = dQ_dVa
        end
    end
    Z_bb = spzeros(nbus, nbus)
    Z_bb = spzeros(nbus, nbus)
    I_gen = spzeros(nbus, ngen)
    for (i,j,v) in zip(generators.bus, collect(1:ngen), ones(ngen))
    I_gen[i,j] = v
    end
    Z_bg = spzeros(nbus, ngen)
    dP_dVm = J[1:nbus, 1:nbus]
    dP_dVa = J[1:nbus, (nbus+1):(2nbus)]
    dQ_dVm = J[(nbus+1):(2nbus), 1:nbus]
    dQ_dVa = J[(nbus+1):(2nbus), (nbus+1):(2nbus)]

    if full == true
        JJ = [ -I_gen    Z_bg    dP_dVm   dP_dVa   I      Z_bb;
                Z_bg    -I_gen   dQ_dVm   dQ_dVa   Z_bb   I    ]
        return JJ
    else
        JJ = Dict()
        JJ[:dP_dVm] = dP_dVm
        JJ[:dP_dVa] = dP_dVa
        JJ[:dQ_dVm] = dQ_dVm
        JJ[:dQ_dVa] = dQ_dVa
        return JJ
    end
end

## -----------------------------------------------------------------------------
## algebraic (complex, vectorized): http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/makeJac.html
## -----------------------------------------------------------------------------
# function PF(xmodel::AbstractArray, Y::AbstractArray, z_idx::Dict)
#     ## setup
#     Vm = xmodel[z_idx[:Vm]]
#     Va = xmodel[z_idx[:Va]]
#     ngen = length([pars[:Pd_R]; pars[:Pd_G]]); @assert(ngen == length([pars[:Pd_R]; pars[:Qd_G]]))
#     nbus = length(pars[:Vm_RGL]); @assert(nbus == length(pars[:Va_RGL]))
#     @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)
#
#     # branch admitances
#     Vtilde = Vm .* exp.(im .* Va)
#     Itilde = Y*Vtilde
#     Stilde = Vtilde .* conj.(Itilde)
#     return Vtilde, Itilde, Stilde
# end
# function PF(y::AbstractArray, data::Dict, full=false)
#     ## setup
#     Y = data[:Y]
#     jac_RGL_idx = data[:jac_RGL_idx]
#     xmodel = data[:xmodel]
#     ngen = length([pars[:Pd_R]; pars[:Pd_G]]); @assert(ngen == length([pars[:Pd_R]; pars[:Qd_G]]))
#     nbus = length(pars[:Vm_RGL]); @assert(nbus == length(pars[:Va_RGL]))
#     @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)
#
#     y_idx = idx_y(vars, idx_type)
#     xmodel[y_idx] = y
#     return PF(xmodel, Y, idx)
# end
# function PFE!(F::AbstractArray, y::AbstractArray, data::Dict)
#     ## setup
#     idx = data[:idx]
#     FF = data[:FF]
#     @assert(length(F) == length(y))
#     ngen = length(idx[:Pg]); @assert(ngen == length(idx[:Qg]))
#     nbus = length(idx[:Vm]); @assert(nbus == length(idx[:Va]))
#     nload = length(idx[:load])
#     Vtilde, Itilde, Stilde = PF(y, data)
#     FF[1:nbus] = real.(Stilde) - data[:Pnet]
#     FF[(nbus+1):(2nbus)] = imag.(Stilde) - data[:Qnet]
#     F = FF[idx[:f]]
#     nothing
# end
#
# function dStilde_dVtilde(Vtilde::AbstractArray, Y::AbstractArray)
#     Itilde = Y * Vtilde
#     diagV = spdiagm(0 => Vtilde)
#     diagVnorm = spdiagm(0 => Vtilde ./ abs.(Vtilde))
#     diagItilde = spdiagm(0 => Itilde)
#     dStilde_dVm = diagV * conj.(Y * diagVnorm) + conj.(diagItilde) * diagVnorm;
#     dStilde_dVa = im * diagV * conj(diagItilde - Y * diagV);
#     return dStilde_dVm, dStilde_dVa
# end
#
# function jac_x(xmodel::AbstractArray, Y::AbstractArray, x_idx::Dict, matchnumerical=true)
#     Vm = xmodel[x_idx[:Vm]]
#     Va = xmodel[x_idx[:Va]]
#     ngen = length(idx[:Pg]); @assert(ngen == length(idx[:Qg]))
#     nbus = length(idx[:Vm]); @assert(nbus == length(idx[:Va]))
#     @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)
#     Vtilde = Vm .* exp.(im .* Va)
#
#     ## jacobian
#     dStilde_dVm, dStilde_dVa = dStilde_dVtilde(Vtilde, Y)
#     dPdVa = real.(dStilde_dVa)
#     dPdVm = real.(dStilde_dVm)
#     dQdVa = imag(dStilde_dVa)
#     dQdVm = imag(dStilde_dVm)
#     Z_bb = spzeros(nbus, nbus)
#     Z_bb = spzeros(nbus, nbus)
#     I_gen = spzeros(nbus, ngen)
#     for (i,j,v) in zip(, collect(1:ngen), ones(ngen))
#         I_gen[i,j] = v
#     end
#     #spdiagm(0 => ones(nbus))[:, 1:ngen]
#     Z_bg = spzeros(nbus, ngen)
#     if matchnumerical
#         return [ -I_gen    Z_bg    dPdVm   dPdVa   I      Z_bb;
#                   Z_bg    -I_gen   dQdVm   dQdVa   Z_bb   I    ]
#     else
#         return [ dPdVa   dPdVm ;
#                  dQdVa   dQdVm ]
#     end
# end