## -----------------------------------------------------------------------------
## solve
## -----------------------------------------------------------------------------
function acopf_solve(opfmodel::JuMP.Model, opfdata::OPFData, warm_point=false)

  #
  # initial point - needed especially for pegase cases
  #
  if warm_point == false
    Pg0,Qg0,Vm0,Va0 = acopf_initialPt_IPOPT(opfdata)
  else
    Pg0 = warm_point[:Pg]
    Qg0 = warm_point[:Qg]
    Vm0 = warm_point[:Vm]
    Va0 = warm_point[:Va]
  end
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

function acopf_solve(M::OPFModel, opfdata::OPFData, warm_point=false); return OPFModel(acopf_solve(M.m, opfdata, warm_point)..., M.kind); end

function acpf_solve(opfmodel::JuMP.Model, opfdata::OPFData, warm_point=false)

  #
  # initial point - assumed to be accounted for in `opfmodel`
  #
  status = :IpoptInit
  status = solve(opfmodel)

  if status != :Optimal
    println("Could not solve the model to optimality.")
  end
  return opfmodel, status
end

function acpf_solve(M::OPFModel, opfdata::OPFData); return OPFModel(acpf_solve(M.m, opfdata)..., M.kind); end

function scacopf_solve(opfmodel::JuMP.Model, opfdata::OPFData, options::Dict, contingencies::Dict=Dict(), warm_point=false, current_rating_bool::Bool=true)
    options[:current_rating] = current_rating_bool
    nbus = length(opfdata.buses); nline = length(opfdata.lines); ngen = length(opfdata.generators)
    R     = opfdata.bus_ref
    G     = filter(x -> x ∉ R, findall(.!isempty.(opfdata.BusGenerators)))
    L     = findall(isempty.(opfdata.BusGenerators))
    not_R = deleteat!(collect(1:nbus), R)

    #
    # initial point - needed especially for pegase cases
    #
    if warm_point == false
        Pg0,Qg0,Vm0,Va0 = acopf_initialPt_IPOPT(opfdata)
    else
        Pg0 = warm_point[:Pg]
        Qg0 = warm_point[:Qg]
        Vm0 = warm_point[:Vm]
        Va0 = warm_point[:Va]
    end
    ## initial dispatch point
    setvalue(getindex(opfmodel, :Pg), Pg0)
    setvalue(getindex(opfmodel, :Qg), Qg0)
    setvalue(getindex(opfmodel, :Vm), Vm0)
    setvalue(getindex(opfmodel, :Va), Va0)

    dp_0, m_0 = get_dispatch_point(opfdata, options)
    for c_id in keys(contingencies)
        for gi in zip([G; R], [first.(opfdata.BusGenerators[G]); first.(opfdata.BusGenerators[R])])
            g = gi[1]
            i = gi[2]
            setvalue(getindex(opfmodel, Symbol("Qg_$(c_id)"))[g], dp_0[:Qg][i])
        end
        for i in L
            setvalue(getindex(opfmodel, Symbol("Vm_$(c_id)"))[i], dp_0[:Vm][i])
        end
        for i in not_R
            setvalue(getindex(opfmodel, Symbol("Va_$(c_id)"))[i], dp_0[:Va][i])
        end
    end

    status = :IpoptInit
    status = solve(opfmodel)

    if status != :Optimal
        println("Could not solve the model to optimality.")
    end
    return opfmodel, status
end

function scacopf_solve(M::OPFModel, opfdata::OPFData, options::Dict, contingencies::Dict=Dict(), warm_point=false); return OPFModel(scacopf_solve(M.m, opfdata, options, contingencies, warm_point)..., M.kind); end

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

function acopf_initialPt_IPOPT_point(opfdata::MPCCases.OPFData)
  Pg, Qg, Vm, Va = acopf_initialPt_IPOPT(opfdata)
  point = Dict()
  point[:Pg] = Pg
  point[:Qg] = Qg
  point[:Vm] = Vm
  point[:Va] = Va
  return point
end

## -----------------------------------------------------------------------------
## reporting
## -----------------------------------------------------------------------------
function acopf_outputAll(opfmodel::JuMP.Model, kind::Symbol, opfdata::MPCCases.OPFData, lossless::Bool=false, current_rating::Bool=false)
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
    elseif kind == :SC
      optvec=zeros(MathProgBase.numvar(opfmodel))
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
      if current_rating == true
        consRhs = zeros(2*nbus+nflowlim)
      else
        consRhs = zeros(2*nbus+2*nflowlim)
      end
    elseif kind ==:S
      if current_rating == true
        consRhs = zeros(4*nbus+nflowlim)
      else
        consRhs = zeros(4*nbus+2*nflowlim)
      end
    elseif kind == :SC
      consRhs = zeros(MathProgBase.numconstr(opfmodel))
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
        # if current_rating
        #   Ys = 1/((lossless ? 0.0 : lines[l].r) + lines[l].x*im);
        #   flowmax/=abs(Ys)^2
        # end

        if ( (consRhs[idx]+flowmax)  >= (1-within/100)^2*flowmax )
          ## NOTE: printing precision for low line limits may be > 100% b/c IPOPT tolerance.
          ## E.g., IPOPT may give rhs > flowmax by 0.00XXX% due to solution precision
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
function acopf_outputAll(M::OPFModel, opfdata::OPFData, options::Dict=DefaultOptions())
  cr   = options[:current_rating]
  loss = options[:lossless]
  return acopf_outputAll(M.m, M.kind, opfdata, loss, cr)
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

function get_opfmodeldata(opfdata::OPFData, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments())
    # parse options
    lossless       = options[:lossless]
    current_rating = options[:current_rating]
    remove_Bshunt  = options[:remove_Bshunt]
    remove_tap     = options[:remove_tap]
    print_level    = options[:print_level]
    feasibility    = options[:feasibility]
    Pg_hi          = adjustments[:Pg_hi]
    Pg_lo          = adjustments[:Pg_lo]
    Qg_hi          = adjustments[:Qg_hi]
    Qg_lo          = adjustments[:Qg_lo]
    Vm_hi          = adjustments[:Vm_hi]
    Vm_lo          = adjustments[:Vm_lo]
    if lossless && !current_rating
        println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
        current_rating = true
    end

    # shortcuts for compactness
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;

    # branch admitances
    YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA, lossless=lossless, remove_Bshunt=remove_Bshunt, remove_tap=remove_tap)

    Y = computeAdmittanceMatrix(opfdata, options)

    opfmodeldata              = Dict()
    opfmodeldata[:lines]      = deepcopy(lines);
    opfmodeldata[:buses]      = deepcopy(buses);
    opfmodeldata[:generators] = deepcopy(generators);
    opfmodeldata[:baseMVA]    = deepcopy(baseMVA);
    opfmodeldata[:BusIdx]     = deepcopy(busIdx);
    opfmodeldata[:FromLines]  = deepcopy(FromLines);
    opfmodeldata[:ToLines]    = deepcopy(ToLines);
    opfmodeldata[:BusGenerators]  = deepcopy(BusGeners);
    opfmodeldata[:YffR]       = deepcopy(YffR);
    opfmodeldata[:YffI]       = deepcopy(YffI);
    opfmodeldata[:YttR]       = deepcopy(YttR);
    opfmodeldata[:YttI]       = deepcopy(YttI);
    opfmodeldata[:YftR]       = deepcopy(YftR);
    opfmodeldata[:YftI]       = deepcopy(YftI);
    opfmodeldata[:YtfR]       = deepcopy(YtfR);
    opfmodeldata[:YtfI]       = deepcopy(YtfI);
    opfmodeldata[:YshR]       = deepcopy(YshR);
    opfmodeldata[:YshI]       = deepcopy(YshI);
    opfmodeldata[:Y]          = deepcopy(Y);
    return opfmodeldata
end

function update_loadings!(opfdata::OPFData, options::Dict,
                          loading::Float64=DefaultLoading(), adj_pf::Float64=DefaultAdjPF())
    Pd_new, Qd_new = get_loadings(opfdata, options, loading, adj_pf)
    opfdata.buses.Pd .= Pd_new
    opfdata.buses.Qd .= Qd_new
    nothing
end

function update_ratings_max!(opfdata::OPFData, options::Dict)
  Y     = computeAdmittanceMatrix(opfdata, options)
  f     = opfdata.lines.from
  t     = opfdata.lines.to
  Y_tf  = [Y[tt,ff] for (tt,ff) in zip(t,f)]
  Y_ft  = [Y[ff,tt] for (tt,ff) in zip(t,f)]
  V_t   = opfdata.buses.Vmax[t]
  V_f   = opfdata.buses.Vmax[f]
  Yabs2 = max.(abs2.(Y_tf), abs2.(Y_ft))

  if options[:current_rating] == true
    max_ratings = (V_t .^2 .+ V_f .^2 .+ 2 .* (V_t .* V_f)) .* Yabs2
  else
    throw("Max power ratings not yet implemented.")
  end
  opfdata.lines.rateA .= max_ratings
end

function get_loadings(opfdata::OPFData, options::Dict,
                      loading::Float64=DefaultLoading(), adj_pf::Float64=DefaultAdjPF())
    ## get total suppliable Pg & Qg
    total_Pg = sum(opfdata.generators.Pmax) * opfdata.baseMVA
    total_Qg = sum(opfdata.generators.Qmax) * opfdata.baseMVA
    total_Sg = total_Pg + im*total_Qg
    PFg      = total_Pg / abs(total_Sg)  # demand power factor
    total_Pd = sum(opfdata.buses.Pd)
    total_Qd = sum(opfdata.buses.Qd)
    total_Sd = total_Pd + im*total_Qd
    PFd      = total_Pd / abs(total_Sd)  # demand power factor
    cosθ     = (adj_pf == 0.0) ? max(PFd, PFg) : adj_pf
    total_Pd_new = loading * total_Pg
    total_Qd_new = total_Pd_new * tan(acos(cosθ))

    ## modify Pd & Qd to fixed percentage of total suppliable Pg & Qg
    Pd_new = (opfdata.buses.Pd ./ total_Pd) .* total_Pd_new
    Qd_new = (opfdata.buses.Qd ./ total_Qd) .* total_Qd_new
    return Pd_new, Qd_new
end

function get_dispatch_point(opfdata::OPFData, options::Dict, adjustments::Dict=DefaultAdjustments())
    """
    solve an OPF problem (cold) to get a dispatch point
    """
    M = acopf_model(opfdata, options, adjustments)
    M = acopf_solve(M, opfdata)
    x_dispatch = MathProgBase.getsolution(M.m.internalModel)
    dispatch_point = Dict()
    dispatch_point[:Pg] = deepcopy(x_dispatch[[x.col for x in getindex(M.m, :Pg)]])
    dispatch_point[:Qg] = deepcopy(x_dispatch[[x.col for x in getindex(M.m, :Qg)]])
    dispatch_point[:Vm] = deepcopy(x_dispatch[[x.col for x in getindex(M.m, :Vm)]])
    dispatch_point[:Va] = deepcopy(x_dispatch[[x.col for x in getindex(M.m, :Va)]])

    return dispatch_point, M
end

function get_operating_point(dispatch_point::Dict, opfdata::OPFData, options::Dict, tol=DefaultFeasTol())
    """
    solve a PF problem from a dispatch point to get an operating point
    """
    ## TODO: RETURN INTERMEDIATE VALUE IF INFEASIBLE
    M = acpf_model(dispatch_point, opfdata, options, tol)
    M = acpf_solve(M, opfdata)
    operating_point = Dict()
    operating_point[:Pg] = deepcopy(dispatch_point[:Pg])
    operating_point[:Qg] = deepcopy(getvalue(M.m[:Qg]))
    operating_point[:Vm] = deepcopy(getvalue(M.m[:Vm]))
    operating_point[:Va] = deepcopy(getvalue(M.m[:Va]))
    return operating_point, M
end

## -----------------------------------------------------------------------------
## helpers: computation
## -----------------------------------------------------------------------------
function get_flowmag2s(VM::Array{Float64,1}, VA::Array{Float64,1}, Y::AbstractArray, opfdata::OPFData, options::Dict, c::Bool=false)
    lines = opfdata.lines; busIdx = opfdata.BusIdx; nline = length(lines)
    if c == true
      nline_orig = nline + 1
    else
      nline_orig = nline
    end
    current2s = zeros(nline_orig)  ## NOTE: HARDCODED 1 since only one asset removed
    for l in 1:nline_orig
        line = lines[(lines.id .== l)]
        if !isempty(line)
            line = first(line)
            f = line.from; t = line.to
            Y_tf = Y[t, f]; Y_ft = Y[f, t]
            ## NOTE: current from Frank & Rebennack OPF primer: eq 5.11 where turns/tap ratios are accounted for in `Y`
            Vm_f = VM[busIdx[f]]; Va_f = VA[busIdx[f]]
            Vm_t = VM[busIdx[t]]; Va_t = VA[busIdx[t]]
            Yabs2 = max(abs2(Y_tf), abs2(Y_ft))
            current2 = Vm_f^2 + Vm_t^2 - 2 * Vm_f * Vm_t * cos(Va_f - Va_t)
            current2 *= Yabs2
            current2s[l] = current2
        end
    end
    return (flowmag2=current2s, id=collect(1:nline_orig))
end

function get_flowmag2s(point::Dict, opfdata::OPFData, options::Dict, c::Bool=false)
    VM = point[:Vm]
    VA = point[:Va]
    Y = computeAdmittanceMatrix(opfdata, options)
    return get_flowmag2s(VM, VA, Y, opfdata, options, c)
end

function get_flowmag2s(M::OPFModel, opfdata::OPFData, options::Dict, c::Bool=false)
    VM = getvalue(M.m[:Vm])
    VA = getvalue(M.m[:Va])
    Y = computeAdmittanceMatrix(opfdata, options)
    return get_flowmag2s(VM, VA, Y, opfdata, options, c)
end

function get_ratings(flowmag2s::Array{Float64,1}, baseMVA::Float64=100.0)
    return sqrt.(flowmag2s) * baseMVA
end

## -----------------------------------------------------------------------------
## helpers: topology
## -----------------------------------------------------------------------------
function get_nonislanding_lines(opfdata::OPFData, options::Dict)
    nonislanding_lines = Int64[]
    for l in eachindex(opfdata.lines)
        opfd = deepcopy(opfdata)
        remove_line!(opfd, l)
        Y = sparse(computeAdmittanceMatrix(opfd, options))
        m = strong_components_map(Y)
        if length(unique(m)) == 1; push!(nonislanding_lines, l); end
    end
    return nonislanding_lines
end

function remove_line!(opfdata::OPFData, l::Int64, verb::Bool=false)
    lines = [x for x in opfdata.lines]
    removed = l ∉ eachindex(opfdata.lines)
    if !removed
        if verb; println("removing line $l"); end
        rl = splice!(lines, l)
        opfdata.lines = StructArray(lines)

        ## adjust FromLines & ToLines
        FromLines, ToLines = mapLinesToBuses(opfdata.buses, opfdata.lines, opfdata.BusIdx)
        opfdata.FromLines .= FromLines
        opfdata.ToLines   .= ToLines
        return rl
    else
        throw("Line $l has already been removed")
    end
end


function reinstate_line!(opfdata::OPFData, l::Int64, rl::MPCCases.Line, verb::Bool=false)
    lines = [x for x in opfdata.lines]
    redundant = any([(rl == x) for x in lines])
    if !redundant
        if verb; println("reinstating line $l"); end
        splice!(lines, (l):(l-1), [rl])
        opfdata.lines = StructArray(lines)

        ## adjust FromLines & ToLines
        FromLines, ToLines = mapLinesToBuses(opfdata.buses, opfdata.lines, opfdata.BusIdx)
        opfdata.FromLines .= FromLines
        opfdata.ToLines   .= ToLines
    else
        throw("Line $l has already been reinstated")
    end
    nothing
end


function set_initial_limits!(opfdata::OPFData, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments())
    ## set ratings at maximum ratings
    update_ratings_max!(opfdata, options)
    max_ratings = deepcopy(opfdata.lines.rateA)

    ## get x̄₀ (base dispatch point)
    dp_0, M_0 = get_dispatch_point(opfdata, options, adjustments)

    ## update ratings to base dispatch point
    flowmag2s_0 = get_flowmag2s(dp_0, opfdata, options)
    ratings_0   = get_ratings(flowmag2s_0, opfdata.baseMVA)
    opfdata.lines.rateA .= ratings_0
    nothing
end

function get_all_contingencies(opfdata::OPFData, options::Dict=DefaultOptions())
    contingencies = Dict()
    nline = length(opfdata.lines)
    for l in 1:nline
        contingencies[l] = (c_type=:line, asset=deepcopy(opfdata.lines[l]))
    end
    return contingencies
end

function update_ratings_flowviol!(point::Dict, opfdata::OPFData, options::Dict,
                                  feas_tol=DefaultFeasTol(), buffer=DefaultBuffer())
    """
    modify `opfdata.lines.rateA` so that `point` is feasibile
    """
    feas, infeas_dict = check_feasibility(point, opfdata, options, feas_tol)
    if !isempty(infeas_dict[:flows])
        ## flow infeasibility
        lines = infeas_dict[:flows]
        flowmag2s   = get_flowmag2s(point, opfdata, options)
        ratings     = get_ratings(flowmag2s, opfdata.baseMVA)
        adj_ratings = max.(ratings, opfdata.lines.rateA)
        opfdata.lines.rateA[lines] .= adj_ratings[lines] .* (1.0 + buffer)
    end
end

function check_feasibility(check_point::Dict, opfdata::OPFData, options::Dict, feas_tol=DefaultFeasTol())
    """
    check feasibility of `point` and return dictionary of violating buses and lines
    """
    ## process
    PG = check_point[:Pg]
    QG = check_point[:Qg]
    VM = check_point[:Vm]
    VA = check_point[:Va]
    Y  = computeAdmittanceMatrix(opfdata, options)
    Pg_hi = opfdata.generators.Pmax
    Pg_lo = opfdata.generators.Pmin
    Qg_hi = opfdata.generators.Qmax
    Qg_lo = opfdata.generators.Qmin
    Vm_hi = opfdata.buses.Vmax
    Vm_lo = opfdata.buses.Vmin
    Va_hi = pi
    Va_lo = -pi
    flowmag2s = get_flowmag2s(VM, VA, Y, opfdata, options)
    flowmax   = (opfdata.lines.rateA ./ opfdata.baseMVA).^2

    ## check componentwise feasibility
    feas  = true
    PG_hi = PG .<= Pg_hi .+ feas_tol
    PG_lo = PG .>= Pg_lo .- feas_tol
    QG_hi = QG .<= Qg_hi .+ feas_tol
    QG_lo = QG .>= Qg_lo .- feas_tol
    VM_hi = VM .<= Vm_hi .+ feas_tol
    VM_lo = VM .>= Vm_lo .- feas_tol
    VA_hi = VA .<= Va_hi .+ feas_tol
    VA_lo = VA .>= Va_lo .- feas_tol
    flows = flowmag2s .<= flowmax .+ feas_tol
    feas *= prod(PG_hi); feas *= prod(PG_lo)
    feas *= prod(QG_hi); feas *= prod(QG_lo)
    feas *= prod(VM_hi); feas *= prod(VM_lo)
    feas *= prod(VA_hi); feas *= prod(VA_lo)
    feas *= prod(flows)

    ## infeasibility results
    infeas_dict = Dict()
    infeas_dict[:PG_hi] = findall(.!PG_hi)
    infeas_dict[:PG_lo] = findall(.!PG_lo)
    infeas_dict[:QG_hi] = findall(.!QG_hi)
    infeas_dict[:QG_lo] = findall(.!QG_lo)
    infeas_dict[:VM_hi] = findall(.!VM_hi)
    infeas_dict[:VM_lo] = findall(.!VM_lo)
    infeas_dict[:VA_hi] = findall(.!VA_hi)
    infeas_dict[:VA_lo] = findall(.!VA_lo)
    infeas_dict[:flows] = findall(.!flows)

    if feas == true
        return true, Dict()
    else
        return false, infeas_dict
    end
end
