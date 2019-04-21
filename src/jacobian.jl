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
    dP_dVm = zeros(nbus, nbus)
    dP_dVa = zeros(nbus, nbus)
    dQ_dVm = zeros(nbus, nbus)
    dQ_dVa = zeros(nbus, nbus)

    ## compute
    for i = 1:nbus # P, Q; equations
        for k = 1:nbus # Vm, Va; buses
            i = busIdx[mod1(i, nbus)]
            k = busIdx[mod1(k, nbus)]
            if i == k
                IDX = [busIdx[x] for x in Y[i,:].nzind]
                P = Vm[i] * sum(Vm[kk] * ( G[i,kk] * cos(Va[i]-Va[kk]) + B[i,kk] * sin(Va[i]-Va[kk]) ) for kk in IDX)
                Q = Vm[i] * sum(Vm[kk] * ( G[i,kk] * sin(Va[i]-Va[kk]) - B[i,kk] * cos(Va[i]-Va[kk]) ) for kk in IDX)
                dP_dVa_ik = -Q       - B[i,i] * Vm[i]^2
                dP_dVm_ik =  P/Vm[i] + G[i,i] * Vm[i]
                dQ_dVa_ik =  P       - G[i,i] * Vm[i]^2
                dQ_dVm_ik =  Q/Vm[i] - B[i,i] * Vm[i]
            else
                dP_dVa_ik =  Vm[i] * Vm[k] * ( G[i,k] * sin(Va[i]-Va[k]) - B[i,k] * cos(Va[i]-Va[k]) )
                dP_dVm_ik =  Vm[i]         * ( G[i,k] * cos(Va[i]-Va[k]) + B[i,k] * sin(Va[i]-Va[k]) )
                dQ_dVa_ik = -Vm[i] * Vm[k] * ( G[i,k] * cos(Va[i]-Va[k]) + B[i,k] * sin(Va[i]-Va[k]) )
                dQ_dVm_ik =  Vm[i]         * ( G[i,k] * sin(Va[i]-Va[k]) - B[i,k] * cos(Va[i]-Va[k]) )
            end
            dP_dVa[i, k] = dP_dVa_ik
            dP_dVm[i, k] = dP_dVm_ik
            dQ_dVa[i, k] = dQ_dVa_ik
            dQ_dVm[i, k] = dQ_dVm_ik
        end
    end
    Z_bb = spzeros(nbus, nbus)
    Z_bb = spzeros(nbus, nbus)
    I_gen = spzeros(nbus, ngen)
    for (i,j,v) in zip(genbus, collect(1:ngen), ones(ngen))
        I_gen[i,j] = v
    end
    Z_bg = spzeros(nbus, ngen)

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

function jac_z_alg(opfmodel_z::AbstractArray,
                   Y::AbstractArray,
                   BusIdx::Dict, BusGenerators::Array{Array{Int64,1},1}, z_idx::Dict, m_idx::Dict,
                   matchnumerical=true)
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
    dP_dPg = zeros(nbus, ngen); dP_dQg = zeros(nbus, ngen)
    dP_dVm = zeros(nbus, nbus)
    dP_dVa = zeros(nbus, nbus)
    dP_dPd = zeros(nbus, nbus); dP_dQd = zeros(nbus, nbus)

    dQ_dPg = zeros(nbus, ngen); dQ_dQg = zeros(nbus, ngen)
    dQ_dVm = zeros(nbus, nbus)
    dQ_dVa = zeros(nbus, nbus)
    dQ_dPd = zeros(nbus, nbus); dQ_dQd = zeros(nbus, nbus)

    for i = 1:nbus # P, Q; equations
        for k = 1:nbus # Vm, Va; buses
            i = BusIdx[mod1(i, nbus)]
            k = BusIdx[mod1(k, nbus)]
            if i == k
                IDX = [BusIdx[x] for x in Y[i,:].nzind]
                P = Vm[i] * sum(Vm[kk] * ( G[i,kk] * cos(Va[i]-Va[kk]) + B[i,kk] * sin(Va[i]-Va[kk]) ) for kk in IDX)
                Q = Vm[i] * sum(Vm[kk] * ( G[i,kk] * sin(Va[i]-Va[kk]) - B[i,kk] * cos(Va[i]-Va[kk]) ) for kk in IDX)
                if !isempty(BusGenerators[i]); dP_dPg[i, k] = -1.0; end
                dP_dPd[i, k] = 1.0
                dP_dVa[i, k] = Y[i,k] == 0 ? 0.0 : -Q       - B[i,i] * Vm[i]^2
                dP_dVm[i, k] = Y[i,k] == 0 ? 0.0 :  P/Vm[i] + G[i,i] * Vm[i]
                dQ_dVa[i, k] = Y[i,k] == 0 ? 0.0 :  P       - G[i,i] * Vm[i]^2
                dQ_dVm[i, k] = Y[i,k] == 0 ? 0.0 :  Q/Vm[i] - B[i,i] * Vm[i]
                if !isempty(BusGenerators[i]); dQ_dQg[i, k] = -1.0; end
                dQ_dQd[i, k] = 1.0
            else
                dP_dVa[i, k] = Y[i,k] == 0 ? 0.0 :  Vm[i] * Vm[k] * ( G[i,k] * sin(Va[i]-Va[k]) - B[i,k] * cos(Va[i]-Va[k]) )
                dP_dVm[i, k] = Y[i,k] == 0 ? 0.0 :  Vm[i]         * ( G[i,k] * cos(Va[i]-Va[k]) + B[i,k] * sin(Va[i]-Va[k]) )
                dQ_dVa[i, k] = Y[i,k] == 0 ? 0.0 : -Vm[i] * Vm[k] * ( G[i,k] * cos(Va[i]-Va[k]) + B[i,k] * sin(Va[i]-Va[k]) )
                dQ_dVm[i, k] = Y[i,k] == 0 ? 0.0 :  Vm[i]         * ( G[i,k] * sin(Va[i]-Va[k]) - B[i,k] * cos(Va[i]-Va[k]) )
            end
        end
    end
    #### aggregate components
    J = [ dP_dPg   dP_dQg   dP_dVm   dP_dVa   dP_dPd   dP_dQd;
          dQ_dPg   dQ_dQg   dQ_dVm   dQ_dVa   dQ_dPd   dQ_dQd ]
    JJ = Dict()
    JJ[:dP_dPg] = dP_dPg
    JJ[:dP_dQg] = dP_dQg
    JJ[:dP_dVm] = dP_dVm
    JJ[:dP_dVa] = dP_dVa
    JJ[:dP_dPd] = dP_dPd
    JJ[:dP_dQd] = dP_dQd
    JJ[:dQ_dPg] = dQ_dPg
    JJ[:dQ_dQg] = dQ_dQg
    JJ[:dQ_dVm] = dQ_dVm
    JJ[:dQ_dVa] = dQ_dVa
    JJ[:dQ_dPd] = dQ_dPd
    JJ[:dQ_dQd] = dQ_dQd
    #### partition aggregated Jacobian
    dF = Dict()
    dF[:dF_dx] = J[m_idx[:F], m_idx[:x]]
    dF[:dF_dy] = J[m_idx[:F], m_idx[:y]]

    ## return
    if matchnumerical == true
        return J
    else
        return J, JJ, dF
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