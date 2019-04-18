## -----------------------------------------------------------------------------
## numerical
## -----------------------------------------------------------------------------
function jac_x(x::Array{Float64,1}; model::JuMP.NLPEvaluator)
    nconstr = MathProgBase.numconstr(model.m)
    nvar = MathProgBase.numvar(model.m)
    J = spzeros(nconstr, nvar)
    j!(J, x, model=model)
    return J
end
function dFdy_dFdx(x::Array{Float64,1}, idx::Dict; model::JuMP.NLPEvaluator)
    ## evaluate
    J = jac_x(x, model=model)

    ## return
    dFdy = J[idx[:f], idx[:y]]
    dFdx = J[idx[:f], idx[:x]]
    return dFdy, dFdx
end

function dFdy_dFdx(xmodel::AbstractArray, Y::AbstractArray, idx::Dict, matchnumerical=true)
    ## evaluate
    J = jac_x(xmodel, Y, idx, matchnumerical)

    ## return
    dFdy = J[idx[:f], idx[:y]]
    dFdx = J[idx[:f], idx[:x]]
    return dFdy, dFdx
end

function get_idx_sets(opfdata::OPFData; full::Bool=true)
    ngen = length(opfdata.generators); nbus = length(opfdata.buses)
    busIdx = opfdata.BusIdx

    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus

    ref_idx = busIdx[opfdata.bus_ref]
    gen_idx = filter(x -> busIdx[x] ∉ ref_idx, [x[1] for x in opfdata.BusGenerators[opfdata.generators.bus]])
    load_idx = filter(x -> x ∉ Set([gen_idx; ref_idx]), [busIdx[x] for x in opfdata.buses.bus_i])  ## purely load buses

    if full == true
        #        Qg^{G ∪ R}                           Vm^L                       Va^{L ∪ G}
        y_idx = [Qg_idx_offset .+ [gen_idx; ref_idx]; Vm_idx_offset .+ load_idx; Va_idx_offset .+ [gen_idx; load_idx]]
        #        Pd^{R ∪ G ∪ L}                    Qd^{R ∪ G ∪ L}
        x_idx = [Pd_idx_offset .+ collect(1:nbus); Qd_idx_offset .+ collect(1:nbus)]
        #        P^{G ∪ L}          Q^{R ∪ G ∪ L}
        f_idx = [gen_idx; load_idx; nbus .+ [ref_idx; gen_idx; load_idx]]
    elseif full == false
        #        Vm^L                       Va^{L ∪ G}
        y_idx = [Vm_idx_offset .+ load_idx; Va_idx_offset .+ [gen_idx; load_idx]]
        #        Pd^{R ∪ G ∪ L}                    Qd^{R ∪ G ∪ L}
        x_idx = [Pd_idx_offset .+ collect(1:nbus); Qd_idx_offset .+ collect(1:nbus)]
        #        P^{G ∪ L}          Q^{R ∪ G ∪ L}
        f_idx = [gen_idx; load_idx; nbus .+ load_idx]
    end

    idx = Dict()
    idx[:ref] = ref_idx
    idx[:gen] = gen_idx
    idx[:load] = load_idx
    idx[:y] = y_idx
    idx[:x] = x_idx
    idx[:f] = f_idx
    idx[:Pg] = Pg_idx_offset .+ collect(1:ngen)
    idx[:Qg] = Qg_idx_offset .+ collect(1:ngen)
    idx[:Vm] = Vm_idx_offset .+ collect(1:nbus)
    idx[:Va] = Va_idx_offset .+ collect(1:nbus)
    idx[:Pd] = Pd_idx_offset .+ collect(1:nbus)
    idx[:Qd] = Qd_idx_offset .+ collect(1:nbus)
    idx[:gen_eqs] = opfdata.generators.bus
    if full == true
        idx_type = :full
    else
        idx_type = :red
    end
    idx[:type] = idx_type
    return idx
end

## -----------------------------------------------------------------------------
## analytic: http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/makeJac.html
## -----------------------------------------------------------------------------
function get_values(M::OPFModel)
  @assert(M.kind == :S)
  values = Dict{Symbol, Array{Float64,1}}()
  values[:Pg] = getvalue(M.m[:Pg])
  values[:Qg] = getvalue(M.m[:Qg])
  values[:Vm] = getvalue(M.m[:Vm])
  values[:Va] = getvalue(M.m[:Va])
  values[:Pd] = getvalue(M.m[:Pd])
  values[:Qd] = getvalue(M.m[:Qd])
  values[:x] = [values[:Pg]; values[:Qg]; values[:Vm]; values[:Va]; values[:Pd]; values[:Qd]]
  return values
end

function PQnet(opfmodel::OPFModel, opfdata::OPFData)
    nbus = length(opfdata.buses)
    Pd = opfdata.buses.Pd ./ 100.0
    Qd = opfdata.buses.Qd ./ 100.0
    Pg = zeros(nbus)
    Qg = zeros(nbus)
    Pg[opfdata.generators.bus] = opfdata.generators.Pg
    Qg[opfdata.generators.bus] = opfdata.generators.Qg
    Pg[opfdata.generators.bus] = getvalue(opfmodel.m[:Pg])
    Qg[opfdata.generators.bus] = getvalue(opfmodel.m[:Qg])
    Pnet = Pg - Pd
    Qnet = Qg - Qd
    return Pnet, Qnet
end

function PF(xmodel::AbstractArray, Y::AbstractArray, idx::Dict)
    # setup
    Vm = xmodel[idx[:Vm]]
    Va = xmodel[idx[:Va]]
    ngen = length(idx[:Pg]); @assert(ngen == length(idx[:Qg]))
    nbus = length(idx[:Vm]); @assert(nbus == length(idx[:Va]))
    @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)

    # branch admitances
    Vtilde = Vm .* exp.(im .* Va)
    Itilde = Y*Vtilde
    Stilde = Vtilde .* conj.(Itilde)
    return Vtilde, Itilde, Stilde
end

function PF(y::AbstractArray, data::Dict)
    Y = data[:Y]
    idx = data[:idx]
    xmodel = data[:xmodel]
    Pnet = data[:Pnet]
    Qnet = data[:Qnet]
    ngen = length(idx[:Pg]); @assert(ngen == length(idx[:Qg]))
    nbus = length(idx[:Vm]); @assert(nbus == length(idx[:Va]))
    @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)
    xmodel[idx[:y]] = y
    return PF(xmodel, Y, idx)
end

function PFE!(F::AbstractArray, y::AbstractArray, data::Dict)
    idx = data[:idx]
    ngen = length(idx[:Pg]); @assert(ngen == length(idx[:Qg]))
    nbus = length(idx[:Vm]); @assert(nbus == length(idx[:Va]))
    Vtilde, Itilde, Stilde = PF(y, data)
    nonref = [idx[:gen]; idx[:load]]
    F[1:(nbus-1)] = real.(Stilde)[] - data[:Pnet]
    F[(nbus+1):(2nbus)] = imag.(Stilde) - data[:Qnet]
    nothing
end

function dStilde_dV(Vtilde::AbstractArray, Y::AbstractArray)
    Itilde = Y * Vtilde
    diagV = spdiagm(0 => Vtilde)
    diagVnorm = spdiagm(0 => Vtilde ./ abs.(Vtilde))
    diagItilde = spdiagm(0 => Itilde)
    dStilde_dVm = diagV * conj.(Y * diagVnorm) + conj.(diagItilde) * diagVnorm;
    dStilde_dVa = im * diagV * conj(diagItilde - Y * diagV);
    return dStilde_dVm, dStilde_dVa
end

function jac_x(xmodel::AbstractArray, Y::AbstractArray, idx::Dict, matchnumerical=true)
    Vm = xmodel[idx[:Vm]]
    Va = xmodel[idx[:Va]]
    ngen = length(idx[:Pg]); @assert(ngen == length(idx[:Qg]))
    nbus = length(idx[:Vm]); @assert(nbus == length(idx[:Va]))
    @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)
    Vtilde = Vm .* exp.(im .* Va)

    ## jacobian
    dStilde_dVm, dStilde_dVa = dStilde_dV(Vtilde, Y)
    dPdVa = real.(dStilde_dVa)
    dPdVm = real.(dStilde_dVm)
    dQdVa = imag(dStilde_dVa)
    dQdVm = imag(dStilde_dVm)
    Z_bb = spzeros(nbus, nbus)
    Z_bb = spzeros(nbus, nbus)
    I_gen = spzeros(nbus, ngen)
    for (i,j,v) in zip(idx[:gen_eqs], collect(1:ngen), ones(ngen))
        I_gen[i,j] = v
    end
    #spdiagm(0 => ones(nbus))[:, 1:ngen]
    Z_bg = spzeros(nbus, ngen)
    if matchnumerical
        return [ -I_gen    Z_bg    dPdVm   dPdVa   I      Z_bb;
                  Z_bg    -I_gen   dQdVm   dQdVa   Z_bb   I    ]
    else
        return [ dPdVa   dPdVm ;
                 dQdVa   dQdVm ]
    end
end





G = real.(Y); g_row, g_col, g_val = findnz(G)
  B = imag.(Y); b_row, b_col, b_val = findnz(B)
  for q = 1:nbus # P, Q; equations
    for b = 1:nbus # Vm, Va; buses
     h = busIdx[mod1(b, nbus)]
     k = busIdx[mod1(q, nbus)]
     if h == k
       IDX = Y[h,:].nzind
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
     J[q, b]           = dP_dVm
     J[q, nbus+b]      = dP_dVa
     J[nbus+q, b]      = dQ_dVm
     J[nbus+q, nbus+b] = dQ_dVa
    end
  end
  Z_bb = spzeros(nbus, nbus)
  Z_bb = spzeros(nbus, nbus)
  I_gen = spzeros(nbus, ngen)
  for (i,j,v) in zip(opfdata.generators.bus, collect(1:ngen), ones(ngen))
    I_gen[i,j] = v
  end
  Z_bg = spzeros(nbus, ngen)
  dPdVm = J[1:nbus, 1:nbus]
  dPdVa = J[1:nbus, (nbus+1):(2nbus)]
  dQdVm = J[(nbus+1):(2nbus), 1:nbus]
  dQdVa = J[(nbus+1):(2nbus), (nbus+1):(2nbus)]

  J_numerical__ = J_numerical[:, (2ngen+1):(end)]
  J_numerical__ = J_numerical__[:, 1:2nbus]

  J_numerical__[1:nbus,1:nbus] - dPdVm'
  J_numerical__[1:nbus,(nbus+1):(2nbus)] - dPdVa'
  J_numerical__[(nbus+1):(2nbus), 1:nbus] - dQdVm'
  J_numerical__[(nbus+1):(2nbus), (nbus+1):(2nbus)] - dQdVa'

  JJ = [ -I_gen    Z_bg    dPdVm   dPdVa   I      Z_bb;
          Z_bg    -I_gen   dQdVm   dQdVa   Z_bb   I    ]
  J_numerical = J_numerical[1:60, :]
  JJ - J_numerical


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
               dQ_dVa =  Vm[h] * Vm[k] * ( G[h,k] * cos(Va[h]-Va[k]) + B[h,k] * sin(Va[h]-Va[k]) )
               dQ_dVm = -Vm[h]         * ( G[h,k] * sin(Va[h]-Va[k]) - B[h,k] * cos(Va[h]-Va[k]) )
             end
             # J[q, b]           = dP_dVa
             # J[q, nbus+b]      = dP_dVm
             # J[nbus+q, b]      = dQ_dVa
             # J[nbus+q, nbus+b] = dQ_dVm
             J[q, b]           = dP_dVm
             J[q, nbus+b]      = dP_dVa
             J[nbus+q, b]      = dQ_dVm
             J[nbus+q, nbus+b] = dQ_dVa
           end
         end