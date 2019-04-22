## -----------------------------------------------------------------------------
## numerical
## -----------------------------------------------------------------------------
"""
## `jac_z_num`: compute jacobian of full model (PFE, line-constraints, ...) at full data point `z := (x, y)` from `JuMP`
### arguments:
    - `z::AbstractArray`: full data `z := (x, y) = (Pg; Qg; Vm; Va; Pd; Qd)` from an `OPFModel`
    - `model::JuMP.NLPEvaluator`: NLP evaluator
### returns
    - `J::Array{Float64,2}`: full Jacobian (`dF_dz`)
    `J = [ dP_dPg   dP_dQg   dP_dVm   dP_dVa   dP_dPd   dP_dQd;
          dQ_dPg   dQ_dQg   dQ_dVm   dQ_dVa   dQ_dPd   dQ_dQd
          ...      ...      ...      ...      ...      ...    ]`
"""
function jac_z_num(z::Array{Float64,1}; model::JuMP.NLPEvaluator)
    nconstr = MathProgBase.numconstr(model.m)
    nvar = MathProgBase.numvar(model.m)
    J = spzeros(nconstr, nvar)
    j!(J, z, model=model)
    return J
end

## -----------------------------------------------------------------------------
## algebraic (real, entrywise): http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/makeJac.html
## -----------------------------------------------------------------------------
"""
## `jac_z_alg_ew`: compute (entrywise) jacobian of full PFE at full data point `z := (x, y)` based on http://schevalier.com/wp-content/uploads/2017/02/Power-Flow-and-Covariance-Matrix.pdf
### arguments:
    - `z::AbstractArray`: full data `z := (x, y) = (Pg; Qg; Vm; Va; Pd; Qd)` from an `OPFModel`
    - `Y::AbstractArray`: complex admittance matrix
    - `BusIdx::Dict`: bus ID to IDX map (from `OPFData`)
    - `BusGenerators::Array{Array{Int64,1},1}`: generator ID to bus IDX map (from `OPFData`)
    - `m_idx::Dict`: model index map
### returns
    - `J::Array{Float64,2}`: full Jacobian (`dF_dz`)
    `J = [ dP_dPg   dP_dQg   dP_dVm   dP_dVa   dP_dPd   dP_dQd;
          dQ_dPg   dQ_dQg   dQ_dVm   dQ_dVa   dQ_dPd   dQ_dQd ]`
    - `JJ::Dict`: dictionary of subcomponents of `dF_dz`
    - `dF::Dict`: dictionary of `dF_dx`, `dF_dy`, `dF_dd`
"""
function jac_z_alg_ew(z::AbstractArray, Y::AbstractArray,
                      BusIdx::Dict, BusGenerators::Array{Array{Int64,1},1}, m_idx::Dict)
    ## setup
    #### dimensions
    ngen = length(m_idx[:Pg]); @assert(ngen == length(m_idx[:Qg]))
    nbus = length(m_idx[:Vm]); @assert(nbus == length(m_idx[:Va]))
    @assert(length(z) == 2ngen + 2nbus + 2nbus)
    #### data
    G = real.(Y)
    B = imag.(Y)
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
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
                P = P_i(Vm, Va, i, Y)
                Q = Q_i(Vm, Va, i, Y)
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
    dF[:dF_dd] = J[m_idx[:F], m_idx[:d]]

    ## return
    return J, JJ, dF
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
function dStilde_dVtilde(Vm::AbstractArray, Va::AbstractArray, Y::AbstractArray)
    return dStilde_dVtilde(Vm .* exp.(im .* Va), Y)
end
function dStilde_dVtilde(z::AbstractArray, Y::AbstractArray, m_idx::Dict)
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
    return dStilde_dVtilde(Vm, Va, Y)
end

"""
## `jac_z_alg_vec`: compute (vectorized) jacobian of full PFE at full data point `z := (x, y)` based on http://schevalier.com/wp-content/uploads/2017/02/Power-Flow-and-Covariance-Matrix.pdf
### arguments:
    - `z::AbstractArray`: full data `z := (x, y) = (Pg; Qg; Vm; Va; Pd; Qd)` from an `OPFModel`
    - `Y::AbstractArray`: complex admittance matrix
    - `genbus::Array{Int64,1}`: array of generator buses
    - `m_idx::Dict`: model index map
### returns
    - `J::Array{Float64,2}`: full Jacobian (`dF_dz`)
    `J = [ dP_dPg   dP_dQg   dP_dVm   dP_dVa   dP_dPd   dP_dQd;
          dQ_dPg   dQ_dQg   dQ_dVm   dQ_dVa   dQ_dPd   dQ_dQd ]`
    - `JJ::Dict`: dictionary of subcomponents of `dF_dz`
    - `dF::Dict`: dictionary of `dF_dx`, `dF_dy`, `dF_dd`
"""
function jac_z_alg_vec(z::AbstractArray, Y::AbstractArray,
                       genbus::AbstractArray, m_idx::Dict)
    ngen = length(m_idx[:Pg]); @assert(ngen == length(m_idx[:Qg]))
    nbus = length(m_idx[:Vm]); @assert(nbus == length(m_idx[:Va]))
    @assert(length(z) == 2ngen + 2nbus + 2nbus)
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
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
    Z_bg = spzeros(nbus, ngen)
    J = [ -I_gen    Z_bg    dP_dVm   dP_dVa   I      Z_bb;
           Z_bg    -I_gen   dQ_dVm   dQ_dVa   Z_bb   I    ]
    JJ = Dict()
    JJ[:dP_dPg] = -I_gen
    JJ[:dP_dQg] = Z_bg
    JJ[:dP_dVm] = dP_dVm
    JJ[:dP_dVa] = dP_dVa
    JJ[:dP_dPd] = I
    JJ[:dP_dQd] = Z_bb
    JJ[:dQ_dPg] = Z_bg
    JJ[:dQ_dQg] = -I_gen
    JJ[:dQ_dVm] = dQ_dVm
    JJ[:dQ_dVa] = dQ_dVa
    JJ[:dQ_dPd] = Z_bb
    JJ[:dQ_dQd] = I
    #### partition aggregated Jacobian
    dF = Dict()
    dF[:dF_dx] = J[m_idx[:F], m_idx[:x]]
    dF[:dF_dy] = J[m_idx[:F], m_idx[:y]]
    dF[:dF_dd] = J[m_idx[:F], m_idx[:d]]

    ## return
    return J, JJ, dF
end

## -----------------------------------------------------------------------------
## wrapper
## -----------------------------------------------------------------------------
function jac_z(z::AbstractArray, data::Dict, jac_type::Symbol=:ew)
    ## extract
    if jac_type == :ew
        @assert(haskey(data, :Y))
        @assert(haskey(data, :opfdata) || (haskey(data, :BusIdx) && haskey(data, :BusGenerators)))
        @assert(haskey(data, :m_idx))
        Y             = data[:Y]
        BusIdx        = haskey(data, :BusIdx)        ? data[:BusIdx]        : data[:opfdata].BusIdx
        BusGenerators = haskey(data, :BusGenerators) ? data[:BusGenerators] : data[:opfdata].BusGenerators
        m_idx         = data[:m_idx]
        return jac_z_alg_ew(z, Y, BusIdx, BusGenerators, m_idx)
    elseif jac_type == :vec
        @assert(haskey(data, :Y))
        @assert(haskey(data, :genbus) || haskey(data, :opfdata))
        @assert(haskey(data, :m_idx))
        Y             = data[:Y]
        genbus        = haskey(data, :genbus) ? data[:genbus] : data[:opfdata].generators.bus
        m_idx         = data[:m_idx]
        jac_z_alg_vec(z, Y, genbus, m_idx)
    elseif jac_type == :num
        @assert(haskey(data, :model))
        model         = data[:model]
        jac_z_num(z, model=model)
    end
end