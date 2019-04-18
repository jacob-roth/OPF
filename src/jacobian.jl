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
    yidx = om_y_RGL_idx(jac_RGL_idx, full)
    xidx = om_x_RGL_idx(jac_RGL_idx)
    fidx = om_f_RGL_idx(jac_RGL_idx, full)
    @assert(length(yidx) == length(fidx))
    dFdy = J[fidx, yidx]
    dFdx = J[fidx, xidx]
    return dFdy, dFdx
end

## -----------------------------------------------------------------------------
## numerical
## -----------------------------------------------------------------------------
function jac_z_num(opfmodel_z::Array{Float64,1}; model::JuMP.NLPEvaluator)
    nconstr = MathProgBase.numconstr(model.m)
    nvar = MathProgBase.numvar(model.m)
    J = spzeros(nconstr, nvar)
    j!(J, opfmodel_z, model=model)
    return J
end

## -----------------------------------------------------------------------------
## algebraic (real, entrywise): http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/makeJac.html
## -----------------------------------------------------------------------------
function jac_z_alg_ew(opfmodel_z::AbstractArray,
                      Y::AbstractArray,
                      busIdx::Dict, z_idx::Dict,
                      genbus::AbstractArray,
                      matchnumerical=true)
    """ based on: http://schevalier.com/wp-content/uploads/2017/02/Power-Flow-and-Covariance-Matrix.pdf """
    ## setup
    #### dimensions
    ngen = length(z_idx[:Pg]); @assert(ngen == length(z_idx[:Qg]))
    nbus = length(z_idx[:Vm]); @assert(nbus == length(z_idx[:Va]))
    @assert(length(opfmodel_z) == 2ngen + 2nbus + 2nbus)
    #### data
    G = real.(Y)
    B = imag.(Y)
    Vm = opfmodel_z[z_idx[:Vm]]
    Va = opfmodel_z[z_idx[:Va]]
    #### structure
    J = zeros(2nbus, 2nbus)

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
    for (i,j,v) in zip(genbus, collect(1:ngen), ones(ngen))
        I_gen[i,j] = v
    end
    Z_bg = spzeros(nbus, ngen)
    dP_dVm = J[1:nbus, 1:nbus]
    dP_dVa = J[1:nbus, (nbus+1):(2nbus)]
    dQ_dVm = J[(nbus+1):(2nbus), 1:nbus]
    dQ_dVa = J[(nbus+1):(2nbus), (nbus+1):(2nbus)]

    if matchnumerical == true
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
function dStilde_dVtilde(Vtilde::AbstractArray, Y::AbstractArray)
    Itilde = Y * Vtilde
    diagV = spdiagm(0 => Vtilde)
    diagVnorm = spdiagm(0 => Vtilde ./ abs.(Vtilde))
    diagItilde = spdiagm(0 => Itilde)
    dStilde_dVm = diagV * conj.(Y * diagVnorm) + conj.(diagItilde) * diagVnorm;
    dStilde_dVa = im * diagV * conj(diagItilde - Y * diagV);
    return dStilde_dVm, dStilde_dVa
end

function jac_z_alg_vec(opfmodel_z::AbstractArray,
                       Y::AbstractArray,
                       z_idx::Dict,
                       genbus::AbstractArray,
                       matchnumerical=true)
    Vm = opfmodel_z[z_idx[:Vm]]
    Va = opfmodel_z[z_idx[:Va]]
    ngen = length(z_idx[:Pg]); @assert(ngen == length(z_idx[:Qg]))
    nbus = length(z_idx[:Vm]); @assert(nbus == length(z_idx[:Va]))
    @assert(length(opfmodel_z) == 2ngen + 2nbus + 2nbus)
    Vtilde = Vm .* exp.(im .* Va)

    ## jacobian
    dStilde_dVm, dStilde_dVa = dStilde_dVtilde(Vtilde, Y)
    dP_dVa = real.(dStilde_dVa)
    dP_dVm = real.(dStilde_dVm)
    dQ_dVa = imag(dStilde_dVa)
    dQ_dVm = imag(dStilde_dVm)
    Z_bb = spzeros(nbus, nbus)
    Z_bb = spzeros(nbus, nbus)
    I_gen = spzeros(nbus, ngen)
    for (i,j,v) in zip(genbus, collect(1:ngen), ones(ngen))
        I_gen[i,j] = v
    end
    #spdiagm(0 => ones(nbus))[:, 1:ngen]
    Z_bg = spzeros(nbus, ngen)
    if matchnumerical == true
        return [ -I_gen    Z_bg    dP_dVm   dP_dVa   I      Z_bb;
                  Z_bg    -I_gen   dQ_dVm   dQ_dVa   Z_bb   I    ]
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
## wrapper
## -----------------------------------------------------------------------------
function jac_z(opfmodel_z::AbstractArray, data::Dict, jac_type::Symbol=:ew, matchnumerical=true)
    ## extract
    if jac_type ∈ [:ew, :vec]
        @assert(haskey(data, :Y))
        @assert(haskey(data, :opfdata))
        @assert(haskey(data, :z_idx))
        Y = data[:Y]
        busIdx = data[:opfdata].BusIdx
        genbus = data[:opfdata].generators.bus
        z_idx = data[:z_idx]
    elseif jac_type == :num
        @assert(haskey(data, :model))
        model = data[:model]
    end

    ## compute
    if jac_type == :ew
        return jac_z_alg_ew(opfmodel_z, Y, busIdx, z_idx, genbus, matchnumerical)
    elseif jac_type == :vec
        return jac_z_alg_vec(opfmodel_z, Y, z_idx, genbus, matchnumerical)
    elseif jac_type == :num
        return jac_z_num(opfmodel_z, model=model)
    end
end