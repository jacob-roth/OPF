## -----------------------------------------------------------------------------
## solve
## -----------------------------------------------------------------------------
function acopf_solve(opfmodel::JuMP.Model, opfdata::OPFData)

  #
  # Initial point - needed especially for pegase cases
  #
  Pg0,Qg0,Vm0,Va0 = acopf_initialPt_IPOPT(opfdata)
  setvalue(getindex(opfmodel, :Pg), Pg0)
  setvalue(getindex(opfmodel, :Qg), Qg0)
  setvalue(getindex(opfmodel, :Vm), Vm0)
  setvalue(getindex(opfmodel, :Va), Va0)
  status = :IpoptInit
  status = solve(opfmodel)

  if status != :Optimal
    println("Could not solve the model to optimality.")
  end
  return opfmodel, status
end
function acopf_solve(M::OPFModel, opfdata::OPFData); return OPFModel(acopf_solve(M.m, opfdata)..., M.kind); end

function ccacopf_solve(opfmodel::JuMP.Model, opfdata::OPFData)

  #
  # Initial point - needed especially for pegase cases
  #
  sm = OPF.sacopf_model(opfdata)
  sm = OPF.acopf_solve(sm, opfdata)
  sm_eval = setup(sm.m);               ## stochastic model evaluator
  sm_zbar = deepcopy(sm_eval.last_x);  ## stochastc model equilibrium z̄
  ## set OPF variables
  Pg_bar = getvalue(getindex(sm.m, :Pg))
  Qg_bar = getvalue(getindex(sm.m, :Qg))
  Vm_bar = getvalue(getindex(sm.m, :Vm))
  Va_bar = getvalue(getindex(sm.m, :Va))
  Pd_bar = getvalue(getindex(sm.m, :Pd))
  Qd_bar = getvalue(getindex(sm.m, :Qd))
  setvalue(getindex(opfmodel, :Pg), Pg_bar)
  setvalue(getindex(opfmodel, :Qg), Qg_bar)
  setvalue(getindex(opfmodel, :Vm), Vm_bar)
  setvalue(getindex(opfmodel, :Va), Va_bar)
  setvalue(getindex(opfmodel, :Pd), Pd_bar)
  setvalue(getindex(opfmodel, :Qd), Qd_bar)
  ## set sensitivitites
  Y = computeAdmittanceMatrix(opfdata)
  m_idx = OPF.model_idx(opfdata)
  z_idx = OPF.om_z_idx(opfdata)
  J, JJ, dF = OPF.jac_z_alg(sm_zbar, Y, opfdata.BusIdx, opfdata.BusGenerators, z_idx, m_idx, false)
  Γ = dF[:dF_dx] \ -dF[:dF_dy]
  ζ = zeros(size(Γ))
  setvalue(getindex(opfmodel, :Gamma), Γ)
  setvalue(getindex(opfmodel, :zeta), ζ)
  println("Setting initial point for CC-ACOPF")
  println("Pg_bar", Pg_bar)
  println("Qg_bar", Qg_bar)
  println("Vm_bar", Vm_bar)
  println("Va_bar", Va_bar)
  println("Pd_bar", Pd_bar)
  println("Qd_bar", Qd_bar)
  status = :IpoptInit
  status = solve(opfmodel)

  if status != :Optimal
    println("Could not solve the model to optimality.")
  end
  return opfmodel, status
end
function ccacopf_solve(M::OPFModel, opfdata::OPFData); return OPFModel(ccacopf_solve(M.m, opfdata)..., M.kind); end

# Compute initial point for IPOPT based on the values provided in the case data
function acopf_initialPt_IPOPT(opfdata::MPCCases.OPFData)
  Pg=zeros(length(opfdata.generators)); Qg=zeros(length(opfdata.generators)); i=1
  for g in opfdata.generators
    # set the power levels in in between the bounds as suggested by matpower
    # (case data also contains initial values in .Pg and .Qg - not used with IPOPT)
    Pg[i]=0.5*(g.Pmax+g.Pmin)
    Qg[i]=0.5*(g.Qmax+g.Qmin)
    i=i+1
  end
  @assert i-1==length(opfdata.generators)

  Vm=zeros(length(opfdata.buses)); i=1;
  for b in opfdata.buses
    # set the ini val for voltage magnitude in between the bounds
    # (case data contains initials values in Vm - not used with IPOPT)
    Vm[i]=0.5*(b.Vmax+b.Vmin);
    i=i+1
  end
  @assert i-1==length(opfdata.buses)

  # set all angles to the angle of the reference bus
  Va = opfdata.buses[opfdata.bus_ref].Va * ones(length(opfdata.buses))

  return Pg,Qg,Vm,Va
end

## -----------------------------------------------------------------------------
## reporting
## -----------------------------------------------------------------------------
function acopf_outputAll(opfmodel::JuMP.Model, kind::Symbol, opfdata::MPCCases.OPFData)
  #shortcuts for compactness
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;

  nbus  = length(buses); nline = length(lines); ngen  = length(generators)

  # OUTPUTING
  println("Objective value: ", getobjectivevalue(opfmodel), "USD/hr")
  VM=getvalue(getindex(opfmodel,:Vm)); VA=getvalue(getindex(opfmodel,:Va))
  PG=getvalue(getindex(opfmodel,:Pg)); QG=getvalue(getindex(opfmodel,:Qg))
  if kind == :S
    PD=getvalue(getindex(opfmodel,:Pd)); QD=getvalue(getindex(opfmodel,:Qd))
  end

  println("============================= BUSES ==================================")
  println("  BUS    Vm     Va    |   Pg (MW)    Qg(MVAr) ")   # |    P (MW)     Q (MVAr)")  #|         (load)   ")

  println("                      |     (generation)      ")
  println("----------------------------------------------------------------------")
  for i in 1:nbus
    @printf("%4d | %6.2f  %6.2f | %s  | \n",
	    buses[i].bus_i, VM[i], VA[i]*180/pi,
	    (length(BusGeners[i])==0) ? "   --          --  " : @sprintf("%7.2f     %7.2f", baseMVA*PG[BusGeners[i][1]], baseMVA*QG[BusGeners[i][1]]))
  end
  println("\n")

  within=20 # percentage close to the limits


  nflowlim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nflowlim += 1
    end
  end

  if nflowlim>0
    println("Number of lines with flow limits: ", nflowlim)
    if kind == :D
      optvec=zeros(2*nbus+2*ngen)
    elseif kind == :S
      optvec=zeros(4*nbus+2*ngen)
    end
    optvec[1:ngen]=PG
    optvec[ngen+1:2*ngen]=QG
    optvec[2*ngen+1:2*ngen+nbus]=VM
    optvec[2*ngen+nbus+1:2*ngen+2*nbus]=VA
    if kind == :S
      optvec[2*ngen+2*nbus+1:2*ngen+3*nbus]=PD
      optvec[2*ngen+3*nbus+1:2*ngen+4*nbus]=QD
    end

    d = JuMP.NLPEvaluator(opfmodel)
    MathProgBase.initialize(d, [:Jac])
    if kind == :D
      consRhs = zeros(2*nbus+2*nflowlim)
    elseif kind ==:S
      consRhs = zeros(4*nbus+2*nflowlim)
    end
    MathProgBase.eval_g(d, consRhs, optvec)
    # d = setup(opfmodel)
    # c!(consRhs, optvec, model=d)

    @printf("================ Lines within %d pct of flow capacity ===================\n", within)
    println("Line   From Bus    To Bus    At capacity")

    nlim=0
    for l in 1:nline
      if lines[l].rateA!=0 && lines[l].rateA<1.0e10
        flowmax=(lines[l].rateA/baseMVA)^2
        idx = 2nbus + 2nlim + 1

        if( (consRhs[idx]+flowmax)  >= (1-within/100)^2*flowmax )
          @printf("%3d      %3d      %3d        %5.3f%%\n", l, lines[l].from, lines[l].to, 100*sqrt((consRhs[idx]+flowmax)/flowmax))
        elseif( (consRhs[idx + 1]+flowmax)  >= (1-within/100)^2*flowmax )
          @printf("%3d      %3d      %3d        %5.3f%%\n", l, lines[l].from, lines[l].to, 100*sqrt((consRhs[idx + 1]+flowmax)/flowmax))
        end
        nlim += 1
      end
    end
  end

  return
end
function acopf_outputAll(M::OPFModel, opfdata::OPFData); return acopf_outputAll(M.m, M.kind, opfdata); end

"""
## `get_values`: get partitioned values of aggregate `OPFModel`'s `z` vector
### arguments:
    - `opfmodel::OPFModel`: opf model
### returns:
    - `values::Dict`: dictionary of partitioned values in `OPFModel`'s order
"""
function get_values(opfmodel::OPFModel)
  @assert(opfmodel.kind == :S)
  values = Dict{Symbol, Array{Float64,1}}()
  values[:Pg] = getvalue(opfmodel.m[:Pg])
  values[:Qg] = getvalue(opfmodel.m[:Qg])
  values[:Vm] = getvalue(opfmodel.m[:Vm])
  values[:Va] = getvalue(opfmodel.m[:Va])
  values[:Pd] = getvalue(opfmodel.m[:Pd])
  values[:Qd] = getvalue(opfmodel.m[:Qd])
  values[:z] = [values[:Pg]; values[:Qg]; values[:Vm]; values[:Va]; values[:Pd]; values[:Qd]]
  return values
end

## -----------------------------------------------------------------------------
## indices
## -----------------------------------------------------------------------------
"""
## `RGL_id`: get ID sets (`bus_i` name sets) for `R`ef, `G`en, and `L`oad buses
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
### returns:
    - `buses_RGL_id::Dict`: dictionary of IDs for `R`, `G`, `L` sets; apply to vector of size `nbus`
    - `gens_RGL_id::Dict`: dictionary of IDs for `R`, `G`, sets; apply to vector of size `ngen`
"""
function RGL_id(opfdata::OPFData)
    #### `opfdata.buses` ordering
    buses_R_id = [opfdata.bus_ref]
    buses_G_id = filter(x -> x ∉ buses_R_id, opfdata.generators.bus)
    buses_L_id = filter(x -> x ∉ Set([buses_G_id; buses_R_id]), opfdata.buses.bus_i)  ## purely load buses
    #### `opfdata.generators` ordering (NOTE: assuming generators' `ID`s are same as index in generator array
    gens_R_id = findall(opfdata.generators.bus .== opfdata.bus_ref)
    gens_G_id = findall(opfdata.generators.bus .!= opfdata.bus_ref)

    buses_RGL_id = Dict()
    buses_RGL_id[:R] = buses_R_id
    buses_RGL_id[:G] = buses_G_id
    buses_RGL_id[:L] = buses_L_id

    gens_RG_id = Dict()
    gens_RG_id[:R] = gens_R_id
    gens_RG_id[:G] = gens_G_id
    return buses_RGL_id, gens_RG_id
end

"""
## `RGL_idx`: get idx sets (index of a bus's `bus_i` ID in `buses`) for `R`ef, `G`en, and `L`oad buses
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
### returns:
    - `b_RGL_idx::Dict`: dictionary of indices for `R`, `G`, `L` sets; apply to vector of size `nbus`
    - `g_RG_idx::Dict`: dictionary of indices for `R`, `G` sets; apply to vector of size `ngen`
"""
function RGL_idx(buses_RGL_id::Dict, gens_RG_id::Dict, busIdx::Dict)
    b_RGL_idx = Dict()
    b_RGL_idx[:R] = [busIdx[x] for x in buses_RGL_id[:R]]
    b_RGL_idx[:G] = [busIdx[x] for x in buses_RGL_id[:G]]
    b_RGL_idx[:L] = [busIdx[x] for x in buses_RGL_id[:L]]

    g_RG_idx = Dict()
    g_RG_idx[:R] = gens_RG_id[:R]
    g_RG_idx[:G] = gens_RG_id[:G]
    return b_RGL_idx, g_RG_idx
end
function RGL_idx(opfdata::OPFData)
    buses_RGL, gens_RG = RGL_id(opfdata)
    return RGL_idx(buses_RGL, gens_RG, opfdata.BusIdx)
end

function model_idx(opfdata::OPFData)
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  BusIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators); nload = length(findall(buses.bustype .== 1))
  opfmodel = Model(solver=IpoptSolver(print_level=0))
  @variable(opfmodel,  generators[i].Pmin  <= Pg[i=1:ngen] <= generators[i].Pmax)
  @variable(opfmodel,  generators[i].Qmin  <= Qg[i=1:ngen] <= generators[i].Qmax)
  @variable(opfmodel,  buses[i].Vmin       <= Vm[i=1:nbus] <= buses[i].Vmax)
  @variable(opfmodel, -pi                  <= Va[i=1:nbus] <= pi)
  @variable(opfmodel,  buses[i].Pd/baseMVA <= Pd[i=1:nbus] <= buses[i].Pd/baseMVA)
  @variable(opfmodel,  buses[i].Qd/baseMVA <= Qd[i=1:nbus] <= buses[i].Qd/baseMVA)
  ## partition variables
  b_RGL_idx, g_RGL_idx = RGL_idx(opfdata)
  #### unknown
  x = [Vm[b_RGL_idx[:L]]; Va[b_RGL_idx[:G]]; Va[b_RGL_idx[:L]]; Qg[g_RGL_idx[:G]]];
  #### control
  u = [Pg[g_RGL_idx[:G]]; Vm[b_RGL_idx[:G]]; Vm[b_RGL_idx[:R]]];
  #### parameter
  p = Va[b_RGL_idx[:R]];
  #### uncertainty
  d = [Pd; Qd];
  #### aggregate "known"
  y = [u; p; d];
  #### dims
  nx = length(x); nu = length(u); np = length(p); nd = length(d); ny = length(y)
  xidx = [xx.col for xx in x]  ## index in model `z`
  uidx = [xx.col for xx in u]  ## index in model `z`
  pidx = [xx.col for xx in p]  ## index in model `z`
  didx = [xx.col for xx in d]  ## index in model `z`
  yidx = [xx.col for xx in y]  ## index in model `z`
  Fidx = [b_RGL_idx[:L]; b_RGL_idx[:G]; nbus .+ b_RGL_idx[:L]; nbus .+ b_RGL_idx[:G]]  ## index in 2nbus equations
  idx = Dict()
  idx[:x] = xidx
  idx[:u] = uidx
  idx[:p] = pidx
  idx[:d] = didx
  idx[:y] = yidx
  idx[:F] = Fidx
  return idx
end







## -----------------------------------------------------------------------------
## to be phased out...
## -----------------------------------------------------------------------------
"""
## `om_z_idx`: get index sets to extract the following values
    - `Pg` and `Qg` in gen-sorted order, i.e., GEN1, GEN2, ... (not BUS sorted order)
    - `Vm`, `Va`, `Pd`, `Qd` in bus-sorted order, i.e., BUS1, BUS2, ...
##  from a vector `om_z` ("opfmodel-z") which contains an `OPFModel`'s aggregate variable denoted by `z`
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
### returns:
    - `z_idx::Dict`: dictionary of index sets for extracting values from `om_z`
### example:
```
julia> z_idx = om_z_idx(opfdata)
julia> om_z = sm_zbar
julia> Vm = om_z[z_idx[:Vm]]
```
"""
function om_z_idx(opfdata::OPFData)
    nbus = length(opfdata.buses)
    ngen = length(opfdata.generators)
    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
    z_idx = Dict()
    z_idx[:Pg] = Pg_idx_offset .+ collect(1:ngen)
    z_idx[:Qg] = Qg_idx_offset .+ collect(1:ngen)
    z_idx[:Vm] = Vm_idx_offset .+ collect(1:nbus)
    z_idx[:Va] = Va_idx_offset .+ collect(1:nbus)
    z_idx[:Pd] = Pd_idx_offset .+ collect(1:nbus)
    z_idx[:Qd] = Qd_idx_offset .+ collect(1:nbus)
    return z_idx
end

"""
## `om_x_RGL_idx`: get `RGL`-partitioned index sets to extract the following
    - `x` (parameters): `Pd_{R,G,L}`, `Qd_{R,G,L}`
##  from a vector `om_z` ("opfmodel-z") which contains an `OPFModel`'s aggregate variable denoted by `z`
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
### returns:
    - `x_RGL_idx::Dict`: dictionary of `RGL`-partitioned index sets for extracting values from `om_z`
"""
function om_x_RGL_idx(opfdata::OPFData)
    ngen = length(opfdata.generators); nbus = length(opfdata.buses)
    busIdx = opfdata.BusIdx

    ## offsets
    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
    bus_idx, gen_idx = RGL_idx(opfdata)

    x_RGL_idx = Dict()
    ## parameters (`x`)
    #### Pd^{R ∪ G ∪ L}
    x_RGL_idx[:Pd_R]   = Pd_idx_offset .+ bus_idx[:R]
    x_RGL_idx[:Pd_G]   = Pd_idx_offset .+ bus_idx[:G]
    x_RGL_idx[:Pd_L]   = Pd_idx_offset .+ bus_idx[:L]
    x_RGL_idx[:Pd_RGL] = [x_RGL_idx[:Pd_R]; x_RGL_idx[:Pd_G]; x_RGL_idx[:Pd_L]]
    ####  Qd^{R ∪ G ∪ L}
    x_RGL_idx[:Qd_R]   = Qd_idx_offset .+ bus_idx[:R]
    x_RGL_idx[:Qd_G]   = Qd_idx_offset .+ bus_idx[:G]
    x_RGL_idx[:Qd_L]   = Qd_idx_offset .+ bus_idx[:L]
    x_RGL_idx[:Qd_RGL] = [x_RGL_idx[:Qd_R]; x_RGL_idx[:Qd_G]; x_RGL_idx[:Qd_L]]
    #### Vm^{R ∪ G}
    x_RGL_idx[:Vm_R]   = Vm_idx_offset .+ bus_idx[:R]
    x_RGL_idx[:Vm_G]   = Vm_idx_offset .+ bus_idx[:G]
    #### Va^{R}
    x_RGL_idx[:Va_R]   = Va_idx_offset .+ bus_idx[:R]
    return x_RGL_idx
end

"""
## `om_y_RGL_idx`: get `RGL`-partitioned index sets to extract the following
    - `y` (variables): `Qg_{R,G}`, `Vm_L`, `Va_{G,L}`
##  from a vector `om_z` ("opfmodel-z") which contains an `OPFModel`'s aggregate variable denoted by `z`
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
### returns:
    - `y_RGL_idx::Dict`: dictionary of `RGL`-partitioned index sets for extracting values from `om_z`
"""
function om_y_RGL_idx(opfdata::OPFData)
    ngen = length(opfdata.generators); nbus = length(opfdata.buses)
    busIdx = opfdata.BusIdx

    ## offsets
    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
    bus_idx, gen_idx = RGL_idx(opfdata)

    y_RGL_idx = Dict()
    ## variables (`y`)
    #### Qg^{R ∪ G}
    y_RGL_idx[:Qg_R]   = Qg_idx_offset .+ gen_idx[:R]
    y_RGL_idx[:Qg_G]   = Qg_idx_offset .+ gen_idx[:G]
    y_RGL_idx[:Qg_RG]  = [y_RGL_idx[:Qg_R]; y_RGL_idx[:Qg_G]]
    #### Vm^L
    y_RGL_idx[:Vm_L]   = Vm_idx_offset .+ bus_idx[:L]
    #### Va^{G ∪ L}
    y_RGL_idx[:Va_G]   = Va_idx_offset .+ bus_idx[:G]
    y_RGL_idx[:Va_L]   = Va_idx_offset .+ bus_idx[:L]
    y_RGL_idx[:Va_GL]  = [y_RGL_idx[:Va_G]; y_RGL_idx[:Va_L]]
    return y_RGL_idx
end

"""
## `om_pfe_RGL_idx`: get `RGL`-partitioned index sets to extract the following PF equations
    - `f`: `P_G` and `P_L`; `Q_R`, `Q_G`, and `Q_L`
##  from an object of dimension `2nbus`
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
### returns:
    - `pfe_RGL_idx::Dict`: dictionary of `RGL`-partitioned index sets for extracting PF equations
"""
function om_pfe_RGL_idx(opfdata::OPFData)
    ngen = length(opfdata.generators); nbus = length(opfdata.buses)
    busIdx = opfdata.BusIdx

    ## offsets
    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
    bus_idx, gen_idx = RGL_idx(opfdata)

    pfe_RGL_idx = Dict()
    ## equations (`f`)
    #### Pnet^{R ∪ G ∪ L}
    pfe_RGL_idx[:P_G] = bus_idx[:G]
    pfe_RGL_idx[:P_L] = bus_idx[:L]
    pfe_RGL_idx[:P_GL] = [pfe_RGL_idx[:P_G]; pfe_RGL_idx[:P_L]]
    #### Qnet^{R ∪ G ∪ L}
    pfe_RGL_idx[:Q_R] = nbus .+ bus_idx[:R]
    pfe_RGL_idx[:Q_G] = nbus .+ bus_idx[:G]
    pfe_RGL_idx[:Q_L] = nbus .+ bus_idx[:L]
    pfe_RGL_idx[:Q_RGL] = [pfe_RGL_idx[:Q_R]; pfe_RGL_idx[:Q_G]; pfe_RGL_idx[:Q_L]]
    return pfe_RGL_idx
end

"""
## `om_jac_RGL_idx`: get `RGL`-partitioned index sets to extract `x`, `y`, and `pfe` indices
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
### returns:
    - `jac_RGL_idx::Dict`: dictionary of `RGL`-partitioned index sets for extracting PF equations
"""
function om_jac_RGL_idx(opfdata::OPFData)
    jac_RGL_idx = Dict()
    jac_RGL_idx[:x] = om_x_RGL_idx(opfdata)
    jac_RGL_idx[:y] = om_y_RGL_idx(opfdata)
    jac_RGL_idx[:f] = om_pfe_RGL_idx(opfdata)
    return jac_RGL_idx
end

"""
## `om_y_RGL_idx`: get `RGL`-partitioned index sets to extract `y` (variables)
### arguments:
    - `jac_RGL_idx::Dict`: full component dictionary from `om_jac_RGL_idx`
    - `full::Bool`: option to return `Qg_RG` as a variable or not
### returns:
    - `yidx::Vector`: array of indixes for `y := ([Qg_RG], Vm_L, Va_GL)`
"""
function om_y_RGL_idx(jac_RGL_idx::Dict, full=false)
    if full
        yidx = [jac_RGL_idx[:y][:Qg_RG]; jac_RGL_idx[:y][:Vm_L]; jac_RGL_idx[:y][:Va_GL]]
    else
        yidx = [jac_RGL_idx[:y][:Vm_L]; jac_RGL_idx[:y][:Va_GL]]
    end
    return yidx
end

"""
## `om_x_RGL_idx`: get `RGL`-partitioned index sets to extract `x` (parameters)
### arguments:
    - `jac_RGL_idx::Dict`: full component dictionary from `om_jac_RGL_idx`
### returns:
    - `xidx::Vector`: array of indixes for `x := (Pd_RGL, Qd_RGL)`
"""
function om_x_RGL_idx(jac_RGL_idx::Dict)
    xidx = [jac_RGL_idx[:x][:Pd_RGL]; jac_RGL_idx[:x][:Qd_RGL]]
    return xidx
end

"""
## `om_f_RGL_idx`: get `RGL`-partitioned index sets to extract `f` (pfe equations)
### arguments:
    - `jac_RGL_idx::Dict`: full component dictionary from `om_jac_RGL_idx`
    - `full::Bool`: option to return equations for `Qg_RG` or not
### returns:
    - `xidx::Vector`: array of indixes for `f := (P_G, P_L, [Q_R, Q_G], Q_L)`
"""
function om_f_RGL_idx(jac_RGL_idx::Dict, full=false)
    if full
        fidx = [jac_RGL_idx[:f][:P_G]; jac_RGL_idx[:f][:P_L]; jac_RGL_idx[:f][:Q_RGL]]
    else
        fidx = [jac_RGL_idx[:f][:P_G]; jac_RGL_idx[:f][:P_L]; jac_RGL_idx[:f][:Q_L]]
    end
    return fidx
end

# """
# ## `opfmodeljac_RGL_idx`: get `RGL`-partitioned index sets to extract the following PF equations
#     - `P_G` and `P_L`
#     - `Q_R`, `Q_G`, and `Q_L`
# ##  from an object of dimension `2nbus`
# ### arguments:
#     - `opfdata::OPFData`: opf data for a particular case
# ### returns:
#     - `ome_RGL_idx::Dict`: dictionary of `RGL`-partitioned index sets for extracting PF equations
# """
# function opfmodeljac_RGL_idx(opfdata::OPFData)
#
# function jac_idx_y(vars::Dict, full::Bool=true)
#     if full == true
#         return [vars[:Qg_RG]; vars[:Vm_L]; vars[:Va_GL]]
#     else
#         return [vars[:Vm_L]; vars[:Va_GL]]
#     end
# end
# function jac_idx_x(pars::Dict)
#     return [pars[:Pd_RGL]; pars[:Qd_RGL]]
# end
# function jac_idx_f(eqns::Dict, full::Bool=true)
#     if full == true
#         return [eqns[:P_GL]; eqns[:Q_RGL]]
#     else
#         return [eqns[:P_GL]; eqns[:Q_L]]
#     end
# end

## -----------------------------------------------------------------------------
## helpers
## -----------------------------------------------------------------------------
"""
## `PQnet`: compute `Pg - Pd` net injections at each bus
### arguments:
    - `opfmodel::OPFModel`: opf model
### returns:
    - `values::Dict`: dictionary of partitioned values in `OPFModel`'s order
"""
# function PQnet(opfmodel::OPFModel, opfdata::OPFData)
function PQnet(opfmodel, opfdata::OPFData)
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

function nonunique(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatedvector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], xs[i])))
            push!(duplicatedvector,xs[i])
        end
    end
    duplicatedvector
end