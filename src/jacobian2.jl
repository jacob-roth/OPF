## -----------------------------------------------------------------------------
## helper
## -----------------------------------------------------------------------------
function get_xmodel_idx(opfdata::OPFData)
    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
    idx = Dict()
    idx[:Pg] = Pg_idx_offset .+ collect(1:ngen)
    idx[:Qg] = Qg_idx_offset .+ collect(1:ngen)
    idx[:Vm] = Vm_idx_offset .+ collect(1:nbus)
    idx[:Va] = Va_idx_offset .+ collect(1:nbus)
    idx[:Pd] = Pd_idx_offset .+ collect(1:nbus)
    idx[:Qd] = Qd_idx_offset .+ collect(1:nbus)
    return idx
end

function RGL_id(opfdata::OPFData)
    """ bus name/ID of `R`ef, `G`en, and `L`oad bus sets """
    ## indexing for `xmodel = [Pg; Qg; Vm; Va; Pd; Qd]
    #### `opfdata.buses` ordering
    buses_R_id = [opfdata.bus_ref]
    buses_G_id = filter(x -> x ∉ buses_R_id, opfdata.generators.bus)
    buses_L_id = filter(x -> x ∉ Set([buses_G_id; buses_R_id]), opfdata.buses.bus_i)  ## purely load buses
    #### `opfdata.generators` ordering (NOTE: assuming generators' `ID`s are same as index in generator array
    gens_R_id = findall(opfdata.generators.bus .== opfdata.bus_ref)
    gens_G_id = findall(opfdata.generators.bus .!= opfdata.bus_ref)

    buses_RGL = Dict()
    buses_RGL[:R] = buses_R_id
    buses_RGL[:G] = buses_G_id
    buses_RGL[:L] = buses_L_id

    gens_RG = Dict()
    gens_RG[:R] = gens_R_id
    gens_RG[:G] = gens_G_id
    return buses_RGL, gens_RG
end
function RGL_idx(buses_RGL::Dict, gens_RG::Dict, busIdx::Dict)
    """ bus index of `R`ef, `G`en, and `L`oad bus sets """
    b_RGL = Dict()
    b_RGL[:R] = [busIdx[x] for x in buses_RGL[:R]]
    b_RGL[:G] = [busIdx[x] for x in buses_RGL[:G]]
    b_RGL[:L] = [busIdx[x] for x in buses_RGL[:L]]

    g_RG = Dict()
    g_RG[:R] = gens_RG[:R]
    g_RG[:G] = gens_RG[:G]
    return b_RGL, g_RG
end
function RGL_idx(opfdata::OPFData)
    buses_RGL, gens_RG = RGL_id(opfdata)
    return RGL_idx(buses_RGL, gens_RG, opfdata.BusIdx)
end

function jac_idx(opfdata::OPFData)
    """ jacobian index for `xmodel` output of `OPFModel` """
    ngen = length(opfdata.generators); nbus = length(opfdata.buses)
    busIdx = opfdata.BusIdx

    ## offsets
    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
    bus_idx, gen_idx = RGL_idx(opfdata)

    ## variables (`y`)
    vars = Dict()
    #### Qg^{R ∪ G}
    vars[:Qg_R] = Qg_idx_offset .+ gen_idx[:R]
    vars[:Qg_G] = Qg_idx_offset .+ gen_idx[:G]
    vars[:Qg_RG] = [vars[:Qg_R]; vars[:Qg_G]]
    #### Vm^L
    vars[:Vm_L] = Vm_idx_offset .+ bus_idx[:L]
    #### Va^{G ∪ L}
    vars[:Va_G] = Va_idx_offset .+ bus_idx[:G]
    vars[:Va_L] = Va_idx_offset .+ bus_idx[:L]
    vars[:Va_GL] = [vars[:Va_G]; vars[:Va_L]]

    ## parameters (`x`)
    pars = Dict()
    #### Pd^{R ∪ G ∪ L}
    pars[:Pd_R] = Pd_idx_offset .+ bus_idx[:R]
    pars[:Pd_G] = Pd_idx_offset .+ bus_idx[:G]
    pars[:Pd_L] = Pd_idx_offset .+ bus_idx[:L]
    pars[:Pd_RGL] = [pars[:Pd_R]; pars[:Pd_G]; pars[:Pd_L]]
    ####  Qd^{R ∪ G ∪ L}
    pars[:Qd_R] = Qd_idx_offset .+ bus_idx[:R]
    pars[:Qd_G] = Qd_idx_offset .+ bus_idx[:G]
    pars[:Qd_L] = Qd_idx_offset .+ bus_idx[:L]
    pars[:Qd_RGL] = [pars[:Qd_R]; pars[:Qd_G]; pars[:Qd_L]]
    #### Vm
    pars[:Vm_R] = Vm_idx_offset .+ bus_idx[:R]
    pars[:Vm_G] = Vm_idx_offset .+ bus_idx[:G]
    pars[:Vm_RGL] = [pars[:Vm_R]; pars[:Vm_G]; vars[:Vm_L]]
    #### Va
    pars[:Va_R] = Va_idx_offset .+ bus_idx[:R]
    pars[:Va_G] = Va_idx_offset .+ bus_idx[:G]
    pars[:Va_RGL] = [pars[:Va_R]; vars[:Va_GL]]

    ## equations (`f`)
    eqns = Dict()
    #### Pnet^{R ∪ G ∪ L}
    eqns[:P_G] = bus_idx[:G]
    eqns[:P_L] = bus_idx[:L]
    eqns[:P_GL] = [eqns[:P_G]; eqns[:P_L]]
    #### Qnet^{R ∪ G ∪ L}
    eqns[:Q_R] = Qd_idx_offset .+ bus_idx[:R]
    eqns[:Q_G] = Qd_idx_offset .+ bus_idx[:G]
    eqns[:Q_L] = Qd_idx_offset .+ bus_idx[:L]
    eqns[:Q_RGL] = [eqns[:Q_R]; eqns[:Q_G]; eqns[:Q_L]]

    return vars, pars, eqns
end

function idx_y(vars::Dict, full::Bool=true)
    if full == true
        return [vars[:Qg_RG]; vars[:Vm_L]; vars[:Va_GL]]
    else
        return [vars[:Vm_L]; vars[:Va_GL]]
    end
end
function idx_x(pars::Dict)
    return [pars[:Pd_RGL]; pars[:Qd_RGL]]
end
function idx_f(eqns::Dict, full::Bool=true)
    if full == true
        return [eqns[:P_GL]; eqns[:Q_RGL]]
    else
        return [eqns[:P_GL]; eqns[:Q_L]]
    end
end

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

## -----------------------------------------------------------------------------
## numerical
## -----------------------------------------------------------------------------
function jac_x(xmodel::Array{Float64,1}; model::JuMP.NLPEvaluator)
    nconstr = MathProgBase.numconstr(model.m)
    nvar = MathProgBase.numvar(model.m)
    J = spzeros(nconstr, nvar)
    j!(J, xmodel, model=model)
    return J
end
function dFdy_dFdx(xmodel::Array{Float64,1}, jac_idx::Dict; model::JuMP.NLPEvaluator)
    ## evaluate
    J = jac_x(xmodel, model=model)

    ## return
    dFdy = J[jac_idx[:f], jac_idx[:y]]
    dFdx = J[jac_idx[:f], jac_idx[:x]]
    return dFdy, dFdx
end
function dFdy_dFdx(xmodel::AbstractArray, Y::AbstractArray, jac_idx::Dict, matchnumerical=true)
    ## evaluate
    J = jac_x(xmodel, Y, idx, matchnumerical)

    ## return
    dFdy = J[jac_idx[:f], jac_idx[:y]]
    dFdx = J[jac_idx[:f], jac_idx[:x]]
    return dFdy, dFdx
end

## -----------------------------------------------------------------------------
## algebraic (complex): http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/makeJac.html
## -----------------------------------------------------------------------------
function PF(xmodel::AbstractArray, Y::AbstractArray, pars::Dict)
    ## setup
    Vm = xmodel[pars[:Vm_RGL]]
    Va = xmodel[pars[:Va_RGL]]
    ngen = length([pars[:Pd_R]; pars[:Pd_G]]); @assert(ngen == length([pars[:Pd_R]; pars[:Qd_G]]))
    nbus = length(pars[:Vm_RGL]); @assert(nbus == length(pars[:Va_RGL]))
    @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)

    # branch admitances
    Vtilde = Vm .* exp.(im .* Va)
    Itilde = Y*Vtilde
    Stilde = Vtilde .* conj.(Itilde)
    return Vtilde, Itilde, Stilde
end
function PF(y::AbstractArray, data::Dict)
    ## setup
    Y = data[:Y]
    idx_type = data[:idx_type]
    vars = data[:vars]
    pars = data[:pars]
    xmodel = data[:xmodel]
    ngen = length([pars[:Pd_R]; pars[:Pd_G]]); @assert(ngen == length([pars[:Pd_R]; pars[:Qd_G]]))
    nbus = length(pars[:Vm_RGL]); @assert(nbus == length(pars[:Va_RGL]))
    @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)

    y_idx = idx_y(vars, idx_type)
    xmodel[y_idx] = y
    return PF(xmodel, Y, idx)
end
function PFE!(F::AbstractArray, y::AbstractArray, data::Dict)
    ## setup
    idx = data[:idx]
    FF = data[:FF]
    @assert(length(F) == length(y))
    ngen = length(idx[:Pg]); @assert(ngen == length(idx[:Qg]))
    nbus = length(idx[:Vm]); @assert(nbus == length(idx[:Va]))
    nload = length(idx[:load])
    Vtilde, Itilde, Stilde = PF(y, data)
    FF[1:nbus] = real.(Stilde) - data[:Pnet]
    FF[(nbus+1):(2nbus)] = imag.(Stilde) - data[:Qnet]
    F = FF[idx[:f]]
    nothing
end

function dStilde_dVtilde(Vtilde::AbstractArray, Y::AbstractArray)
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
    dStilde_dVm, dStilde_dVa = dStilde_dVtilde(Vtilde, Y)
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

## -----------------------------------------------------------------------------
## algebraic (real): http://www.pserc.cornell.edu/matpower/docs/ref/matpower5.0/makeJac.html
## -----------------------------------------------------------------------------
function dStilde_dVtilde(xmodel::AbstractArray, S_i::Int64, V_j::Int64,
                         opfdata::OPFData, admittance::Dict, idx::Dict)
    ## setup
    YffR = admittance[:YffR]; YffI = admittance[:YffI]
    YttR = admittance[:YttR]; YttI = admittance[:YttI]
    YftR = admittance[:YftR]; YftI = admittance[:YftI]
    YtfR = admittance[:YtfR]; YtfI = admittance[:YtfI]
    YshR = admittance[:YshR]; YshI = admittance[:YshI]
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen = length(generators); nload = length(findall(buses.bustype .== 1))
    Vm = xmodel[idx[:Vm]]
    Va = xmodel[idx[:Va]]
    Pg = xmodel[idx[:Pg]]
    Qg = xmodel[idx[:Qg]]
    Pd = xmodel[idx[:Pd]]
    Qd = xmodel[idx[:Qd]]


    ## jacobian (h,k) entry
    h = S_i  ## power index
    k = V_j  ## bus index
    P = ( + ( reduce(+, YffR[l] for l in FromLines[h]; init=0) + reduce(+, YttR[l] for l in ToLines[h]; init=0) + YshR[h] ) * Vm[h]^2
          +   reduce(+, Vm[h]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[h]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[h]-Va[busIdx[lines[l].to]]  )) for l in FromLines[h]; init=0 )
          +   reduce(+, Vm[h]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[h]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[h]-Va[busIdx[lines[l].from]])) for l in ToLines[h]  ; init=0 )
          - ( reduce(+, baseMVA*Pg[g] for g in BusGeners[h]; init=0) - reduce(+, baseMVA*Pd[l] for l in busIdx[h]; init=0) ) / baseMVA)
    Q = ( + ( reduce(+, -YffI[l] for l in FromLines[h]; init=0) + reduce(+, -YttI[l] for l in ToLines[h]; init=0) - YshI[h] ) * Vm[h]^2
          +   reduce(+, Vm[h]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[h]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[h]-Va[busIdx[lines[l].to]]  )) for l in FromLines[h]; init=0 )
          +   reduce(+, Vm[h]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[h]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[h]-Va[busIdx[lines[l].from]])) for l in ToLines[h]  ; init=0 )
          - ( reduce(+, baseMVA*Qg[g] for g in BusGeners[h]; init=0) - reduce(+, baseMVA*Qd[l] for l in busIdx[h]; init=0) ) / baseMVA)
    if h == k
        dP_dVa =
    else
        dP_dVa = Vm[h] * Vm[k] *
    end

    B[i,j] += YftI[l]
    B[j,i] += YtfI[l]
    B[i,i] += YffI[l]
    B[j,j] += YttI[l]

    G[i,j] += YftR[l]
    G[j,i] += YtfR[l]
    G[i,i] += YffR[l]
    G[j,j] += YttR[l]









function PF_real(xmodel::AbstractArray, G_RGL::AbstractArray, B_RGL::AbstractArray, pars::Dict)
    ## setup
    Vm = xmodel[pars[:Vm_RGL]]
    Va = xmodel[pars[:Va_RGL]]
    G = G_RGL; B = B_RGL
    ngen = length([pars[:Pd_R]; pars[:Pd_G]]); @assert(ngen == length([pars[:Pd_R]; pars[:Qd_G]]))
    nbus = length(pars[:Vm_RGL]); @assert(nbus == length(pars[:Va_RGL]))
    @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)

    # branch admitances
    V_R = Vm .* cos.(Va)
    V_I = Vm .* sin.(Va)
    I_R = G*V_R - B*V_I
    I_I = G*V_I + B*V_R
    S_R = V_R .* I_R + V_I .* I_I
    S_I = V_I .* I_R - V_R .* I_I
    Vtilde = (R=V_R, I=V_I)
    Itilde = (R=I_R, I=I_I)
    Stilde = (R=S_R, I=S_I)
    return Vtilde, Itilde, Stilde
end
# function PF_real(y::AbstractArray, data::Dict)
#     ## setup
#     G = data[:G_RGL]
#     B = data[:B_RGL]
#     vars = data[:vars]
#     pars = data[:pars]
#     xmodel = data[:xmodel]
#     ngen = length([pars[:Pd_R]; pars[:Pd_G]]); @assert(ngen == length([pars[:Pd_R]; pars[:Qd_G]]))
#     nbus = length(pars[:Vm_RGL]); @assert(nbus == length(pars[:Va_RGL]))
#     @assert(length(xmodel) == 2ngen + 2nbus + 2nbus)
#
#     idx_y =
#     xmodel[idx[:y]] = y
#     return PF(xmodel, Y, idx)
# end

function dStilde_dVtilde_real(xmodel::AbstractArray, G_RGL::AbstractArray, B_RGL::AbstractArray, pars::Dict)
    Vm = xmodel[pars[:Vm_RGL]]
    Vtilde, Itilde, Stilde = PF_real(xmodel, G_RGL, B_RGL, pars)
    diagV_R = spdiagm(0 => Vtilde.R)
    diagV_I = spdiagm(0 => Vtilde.I)
    diagVnorm_R = spdiagm(0 => Vtilde.R ./ Vm)
    diagVnorm_I = spdiagm(0 => Vtilde.I ./ Vm)
    diagItilde_R = spdiagm(0 => Itilde.R)
    diagItilde_I = spdiagm(0 => Itilde.I)
    GVR_BVI_norm = G * diagVnorm_R - B * diagVnorm_I
    GVI_BVR_norm = G * diagVnorm_I + B * diagVnorm_R
    GVR_BVI = G * diagV_R - B * diagV_I
    GVI_BVR = G * diagV_I + B * diagV_R

    dStilde_dVm_R = diagV_R * GVR_BVI_norm + diagV_I * GVI_BVR_norm + diagItilde_R * diagVnorm_R - diagItilde_I * diagVnorm_I
    dStilde_dVm_I = -diagV_R * GVI_BVR_norm + diagV_I * GVR_BVI_norm + diagItilde_R * diagVnorm_I + diagItilde_I * diagVnorm_R
    dStilde_dVa_R = diagV_R * (diagItilde_I - GVI_BVR) - diagV_I * (diagItilde_R - GVR_BVI)
    dStilde_dVa_I = diagV_R * (diagItilde_R - GVR_BVI) + diagV_I * (diagItilde_I - GVI_BVR)

    dStilde_dVm = (R=dStilde_dVm_R, I=dStilde_dVm_I)
    dStilde_dVa = (R=dStilde_dVa_R, I=dStilde_dVa_I)
    return dStilde_dVm, dStilde_dVa
end