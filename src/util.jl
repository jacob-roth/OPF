## -----------------------------------------------------------------------------
## solve
## -----------------------------------------------------------------------------
function acopf_solve(opfmodel::JuMP.Model, opfdata::OPFData)

  #
  # initial point - needed especially for pegase cases
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
function acopf_outputAll(opfmodel::JuMP.Model, kind::Symbol, opfdata::MPCCases.OPFData, lossless::Bool=false, current_rating::Bool=true)
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
      if options[:current_rating] == true
        consRhs = zeros(2*nbus+nflowlim)
      else
        consRhs = zeros(2*nbus+2*nflowlim)
      end
    elseif kind ==:S
      if options[:current_rating] == true
        consRhs = zeros(4*nbus+nflowlim)
      else
        consRhs = zeros(4*nbus+2*nflowlim)
      end
    end
    MathProgBase.eval_g(d, consRhs, optvec)
    # d = setup(opfmodel)
    # c!(consRhs, optvec, model=d)

    @printf("================ Lines within %d %% of flow capacity ===================\n", within)
    println("Line   From Bus    To Bus    At capacity")

    idx = 2nbus + 1
    for l in 1:nline
      if lines[l].rateA!=0 && lines[l].rateA<1.0e10
        flowmax=(lines[l].rateA/baseMVA)^2
        if current_rating
          Ys = 1/((lossless ? 0.0 : lines[l].r) + lines[l].x*im);
          flowmax/=abs(Ys)^2
        end

        if ( (consRhs[idx]+flowmax)  >= (1-within/100)^2*flowmax )
          @printf("%3d      %3d      %3d        %5.3f%%\n", l, lines[l].from, lines[l].to, 100*sqrt((consRhs[idx]+flowmax)/flowmax))
        elseif !current_rating
          if ( (consRhs[idx + 1]+flowmax)  >= (1-within/100)^2*flowmax )
            @printf("%3d      %3d      %3d        %5.3f%%\n", l, lines[l].from, lines[l].to, 100*sqrt((consRhs[idx + 1]+flowmax)/flowmax))
          end
        end
        idx += current_rating ? 1 : 2
      end
    end
  end

  return
end
function acopf_outputAll(M::OPFModel, opfdata::OPFData, options::Dict); return acopf_outputAll(M.m, M.kind, opfdata,
                                                                               options[:lossless], options[:current_rating]); end

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
    #### `opfdata.generators` ordering (!! NOTE: assuming generators' `ID`s are same as index in generator array !!)
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

"""
## `model_idx`: get idx sets of `OPFModel`
### arguments:
    - `opfdata::OPFData`: opf data for a particular case
    - `xtilde::Bool`: binary option wheere `true` gets full `x̃` (includes Qg), `false` gets reduced `x` (exlcudes Qg)
### returns:
    - `idx::Dict`: dictionary of indices
"""
function model_idx(opfdata::OPFData, xtilde::Bool=false)
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
  if xtilde == true
    x = [Vm[b_RGL_idx[:L]]; Va[b_RGL_idx[:G]]; Va[b_RGL_idx[:L]]; Qg[g_RGL_idx[:G]]];
  else
    x = [Vm[b_RGL_idx[:L]]; Va[b_RGL_idx[:G]]; Va[b_RGL_idx[:L]]];
  end
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
  if xtilde == true
    Fidx = [b_RGL_idx[:L]; b_RGL_idx[:G]; nbus .+ b_RGL_idx[:L]; nbus .+ b_RGL_idx[:G]]  ## index in 2nbus equations
  else
    Fidx = [b_RGL_idx[:L]; b_RGL_idx[:G]; nbus .+ b_RGL_idx[:L]]  ## index in 2nbus equations
  end
  idx = Dict()
  idx[:x] = xidx
  idx[:u] = uidx
  idx[:p] = pidx
  idx[:d] = didx
  idx[:y] = yidx
  idx[:F] = Fidx
  idx[:Pg] = [x.col for x in getindex(opfmodel, :Pg)]
  idx[:Qg] = [x.col for x in getindex(opfmodel, :Qg)]
  idx[:Vm] = [x.col for x in getindex(opfmodel, :Vm)]
  idx[:Va] = [x.col for x in getindex(opfmodel, :Va)]
  idx[:Pd] = [x.col for x in getindex(opfmodel, :Pd)]
  idx[:Qd] = [x.col for x in getindex(opfmodel, :Qd)]
  return idx
end

## -----------------------------------------------------------------------------
## helpers
## -----------------------------------------------------------------------------
"""
## `PQnet`: compute `Pg - Pd` net injections at each bus at an OPF solution (and thus requires `opfmodel`)
### arguments:
    - `opfmodel::OPFModel`: opf model
    - `opfdata::OPFData`: opf data
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
    # Pg[opfdata.generators.bus] = opfdata.generators.Pg
    # Qg[opfdata.generators.bus] = opfdata.generators.Qg
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
