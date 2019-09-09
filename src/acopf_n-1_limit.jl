## -----------------------------------------------------------------------------
## default options
## -----------------------------------------------------------------------------
function DefaultLoading();    return 0.50; end  # 0.85
function DefaultAdjPF();      return 0.0;  end
function DefaultFeasTol();    return 1e-6; end
function DefaultBuffer();     return 0.25; end
function DefaultMaxIter();    return 10;   end
function DefaultMaxLimit();   return 1e7;  end
function DefaultVPct();       return 1.0;  end
function DefaultTightenAmt(); return 1e-3; end
loading  = DefaultLoading()
adj_pf   = DefaultAdjPF()
feas_tol = DefaultFeasTol()
buffer   = DefaultBuffer()
max_iter = DefaultMaxIter()
pct      = DefaultVPct()

## -----------------------------------------------------------------------------
## set ACOPF N-1 limits
## -----------------------------------------------------------------------------
function set_acopf_n1_limits!(opfdata::OPFData, options::Dict,
                              loading=DefaultLoading(), adj_pf=DefaultAdjPF(),
                              feas_tol=DefaultFeasTol(), buffer=DefaultBuffer(),
                              pct=DefaultVPct(), max_iter=DefaultMaxIter())
    #
    # process
    #
    ## shortcuts
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA;
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen = length(generators);

    ## set new loadings
    update_loadings!(opfdata, options, loading, adj_pf)

    ## set all line limits to zero
    opfdata.lines.rateA .= zeros(nline)

    ## get initial IPOPT point (Va = slack, Vm ≈ 1, Pg & Qg from MPC casefile)
    init_point = acopf_initialPt_IPOPT_point(opfdata)

    #
    # N-0
    #
    ## get initial dispatch point and set initial ratings for 0 line outages
    dispatch_point_0, M_0 = get_dispatch_point(opfdata, options); @assert(M_0.status == :Optimal)
    flowmag2s_0 = get_flowmag2s(dispatch_point_0, opfdata, options)
    ratings_0   = get_ratings(flowmag2s_0, opfdata.baseMVA)
    opfdata.lines.rateA .= max(ratings_0, opfdata.lines.rateA)

    ## ensure that cold start under initial dispatch is solvable from initial point
    ensure_cold_solvability!(init_point, opfdata, options, buffer, pct, max_iter)

    length(opfdata.lines)
    opfdata.lines[1]
    opflines = [x for x in opfdata.lines];
    opflines[1]

    ## get non-islanding lines (contingencies)
    nonislanding_lines = get_nonislanding_lines(opfdata, options)
    nonislanding_lines = nonislanding_lines[1:2]

    #
    # check N-1 (for non-islanding lines)
    #
    n1_feas = false
    dispatch_point = dispatch_point_0
    while !n1_feas
        for l in nonislanding_lines
            println("==========     REMOVED LINE $l    ==========")
            rl = remove_line!(opfdata, l)
            operating_point, M = get_operating_point(dispatch_point, opfdata, options)
            feas, infeas_dict = check_feasibility(operating_point, opfdata, options, 1e-3)
            if feas == false
                update_ratings_flowviol!(     operating_point, opfdata, options, feas_tol, buffer)
                update_ratings_boundviol!(    operating_point, opfdata, options, feas_tol, buffer)
                update_adjustments_boundviol!(operating_point, opfdata, options, adjustments, feas_tol)
            end
            reinstate_line!(opfdata, l, rl)
            println("==========    REINSTATED LINE $l  ==========")
        end
        tighten_adjustments!(adjustments)
        dispatch_point, M = get_dispatch_point(opfdata, options, adjustments); @assert(M_0.status == :Optimal)
        n1_feas = feas
    end
end

## -----------------------------------------------------------------------------
## helpers: adjusting
## -----------------------------------------------------------------------------
function ensure_cold_solvability!(init_point::Dict,
                                  opfdata::OPFData, options::Dict,
                                  feas_tol=DefaultFeasTol(), buffer=DefaultBuffer(),
                                  pct=DefaultVPct(), max_iter=DefaultMaxIter())
    ## check initial solvability
    solved, M, calc_point = check_acopf_solvability(init_point, opfdata, options)
    iter = 1
    while solved == false && iter <= max_iter
        println("~~~~~~~~~~ inner solve iteration $iter ~~~~~~~~~~")
        ## get infeasibile point from within solver
        infeas_point = calc_point

        ## adjust ratings so that infeasible point is now feasible
        ensure_ratings_feasibility!(infeas_point, opfdata, options, feas_tol, buffer)

        ## monitor bus voltage magnitudes to adjust line limits
        bus_adj_hi = findall((opfdata.buses.Vmax .- infeas_point[:Vm])  ./ opfdata.buses.Vmax .< pct/100)  ## too close to upper limit
        bus_adj_lo = findall((infeas_point[:Vm]  .- opfdata.buses.Vmin) ./ opfdata.buses.Vmin .< pct/100)  ## too close to lower limit
        line_adj_hi = unique([findall(opfdata.lines.to .∈ Ref(bus_adj_hi));
                              findall(opfdata.lines.from .∈ Ref(bus_adj_hi))])
        line_adj_lo = unique([findall(opfdata.lines.to .∈ Ref(bus_adj_lo));
                              findall(opfdata.lines.from .∈ Ref(bus_adj_lo))])
        opfdata.lines.rateA[line_adj_lo] .*= (1.0 + buffer)
        opfdata.lines.rateA[line_adj_hi] .*= (1.0 + buffer)

        ## solve and update
        solved, M, calc_point = check_acopf_solvability(init_point, opfdata, options)
        iter += 1
        if iter == max_iter
            throw(ErrorException("Maximum number of line rating scalings reached."))
            # @warn "Maximum number of uniform line rating scalings reached."
        end
        if any(opfdata.lines.rateA .> DefaultMaxLimit())
            throw(ErrorException("Maximum line limit increase reached."))
        end
    end
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

function update_adjustments_boundviol!(point::Dict, opfdata::OPFData, options::Dict,
                                       adjustments::Dict,
                                       feas_tol=DefaultFeasTol())
    ## get infeasibility
    feas, infeas_dict = check_feasibility(point, opfdata, options, feas_tol)

    ## voltage bound infeasibility
    [push!(adjustments[:Vm_hi][:i], i) for i in infeas_dict[:VM_hi]];  unique!(adjustments[:Vm_hi][:i])
    [push!(adjustments[:Vm_lo][:i], i) for i in infeas_dict[:VM_lo]];  unique!(adjustments[:Vm_lo][:i])
    # [push!(adjustments[:Va_hi][:i], i) for i in infeas_dict[:VA_hi]];  unique!(adjustments[:Va_hi][:i])
    # [push!(adjustments[:Va_lo][:i], i) for i in infeas_dict[:VA_lo]];  unique!(adjustments[:Va_lo][:i])

    ## generator bound infeasibility
    [push!(adjustments[:Pg_hi][:i], i) for i in infeas_dict[:PG_hi]];  unique!(adjustments[:Pg_hi][:i])
    [push!(adjustments[:Pg_lo][:i], i) for i in infeas_dict[:PG_lo]];  unique!(adjustments[:Pg_lo][:i])
    [push!(adjustments[:Qg_hi][:i], i) for i in infeas_dict[:QG_hi]];  unique!(adjustments[:Qg_hi][:i])
    [push!(adjustments[:Qg_lo][:i], i) for i in infeas_dict[:QG_lo]];  unique!(adjustments[:Qg_lo][:i])
    nothing
end

function update_ratings_boundviol!(point::Dict, opfdata::OPFData, options::Dict,
                                   feas_tol=DefaultFeasTol(), buffer=DefaultBuffer())
    ## get infeasibility
    feas, infeas_dict = check_feasibility(point, opfdata, options, feas_tol)

    ## voltage bound infeasibility
    buses = unique([infeas_dict[:VM_hi];
                    infeas_dict[:VM_lo];
                    infeas_dict[:VA_hi];
                    infeas_dict[:VA_lo];])
    lines = unique([findall(opfdata.lines.to   .∈ Ref(buses));
                    findall(opfdata.lines.from .∈ Ref(buses))])
    opfdata.lines.rateA[lines] .*= (1.0 + buffer)

    ## generator bound infeasibility
    buses = unique([infeas_dict[:PG_hi];
                    infeas_dict[:PG_lo];
                    infeas_dict[:QG_hi];
                    infeas_dict[:QG_lo];])
    lines = unique([findall(opfdata.lines.to   .∈ Ref(buses));
                    findall(opfdata.lines.from .∈ Ref(buses))])
    opfdata.lines.rateA[lines] .*= (1.0 + buffer)
    nothing
end

function tighten_adjustments!(adjustments::Dict, tighten_amt::Float64=DefaultTightenAmt())
    for k in keys(adjustments)
        if !isempty(adjustments[k][:i])
            adjustments[k][:v] += tighten_amt
        end
    end
    nothing
end

# function ensure_contingency_feasibility!(dispatch_point::Dict,
#                                          opfdata::OPFData, options::Dict,
#                                          feas_tol=DefaultFeasTol(), buffer=DefaultBuffer())
#     """
#     for a dispatch point `dispatch_point`, find an operating point given contingency `contingency`,
#     and ensure that the operating point is line-limit feasible by modifying line limits
#     """
#     ## get operating point
#     operating_point, M = get_operating_point(dispatch_point, opfdata, options)
#     cons    = zeros(MathProgBase.numconstr(M.m))
#     pf_eval = setup(M.m)
#     x       = M.m.internalModel.inner.x
#     c!(cons, x, model=pf_eval)
#
#     ## ensure feasibility
#     ensure_ratings_feasibility!(operating_point, opfdata, options, feas_tol, buffer)
# end

function update_loadings!(opfdata::OPFData, options::Dict,
                       loading::Float64=DefaultLoading(), adj_pf::Float64=DefaultAdjPF())
    Pd_new, Qd_new = get_loadings(opfdata, options, loading, adj_pf)
    println("==========     MODIFYING P & Q     ==========")
    opfdata.buses.Pd .= Pd_new
    opfdata.buses.Qd .= Qd_new
    nothing
end

## -----------------------------------------------------------------------------
## helpers: checking
## -----------------------------------------------------------------------------
function check_acopf_solvability(warm_point::Dict, opfdata::OPFData, options::Dict)
    """
    check the solvability of an ACOPF problem initialized at `start_point` and return the calculated point
    """
    ## build acopf model, solve, and get status
    M   = acopf_model(opfdata, options)
    M   = acopf_solve(M, opfdata, warm_point)
    opt = (M.status == :Optimal)

    ## get computed point
    x_calc     = MathProgBase.getsolution(M.m.internalModel)
    calc_point = Dict()
    calc_point[:Pg] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Pg)]])
    calc_point[:Qg] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Qg)]])
    calc_point[:Vm] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Vm)]])
    calc_point[:Va] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Va)]])
    return opt, M, calc_point
end

function check_feasibility(check_point::Dict, opfdata::OPFData, options::Dict, feas_tol=DefaultFeasTol())
    """
    cehck feasibility of `point` and return dictionary of violating buses and lines
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

## -----------------------------------------------------------------------------
## helpers: dispatch (OPF) and operating (PF) points
## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
## helpers: computation
## -----------------------------------------------------------------------------
function get_flowmag2s(VM::Array{Float64,1}, VA::Array{Float64,1}, Y::AbstractArray, opfdata::OPFData, options::Dict)
    lines = opfdata.lines; busIdx = opfdata.BusIdx; nline = length(lines)
    current2s = zeros(nline)
    for l in 1:nline
        line = lines[l]; f = line.from; t = line.to
        Y_tf = Y[t, f]; Y_ft = Y[f, t]
        ## NOTE: current from Frank & Rebennack OPF primer: eq 5.11 where turns/tap ratios are accounted for in `Y`
        Vm_f = VM[busIdx[f]]; Va_f = VA[busIdx[f]]
        Vm_t = VM[busIdx[t]]; Va_t = VA[busIdx[t]]
        Yabs2 = max(abs2(Y_tf), abs2(Y_ft))
        current2 = Vm_f^2 + Vm_t^2 - 2 * Vm_f * Vm_t * cos(Va_f - Va_t)
        current2 *= Yabs2
        current2s[l] = current2
    end
    return current2s
end

function get_flowmag2s(point::Dict, opfdata::OPFData, options::Dict)
    VM = point[:Vm]
    VA = point[:Va]
    Y = computeAdmittanceMatrix(opfdata, options)
    return get_flowmag2s(VM, VA, Y, opfdata, options)
end

function get_flowmag2s(M::OPFModel, opfdata::OPFData, options::Dict)
    VM = getvalue(M.m[:Vm])
    VA = getvalue(M.m[:Va])
    Y = computeAdmittanceMatrix(opfdata, options)
    return get_flowmag2s(VM, VA, Y, opfdata, options)
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
    removed = l ∉ [x.id for x in lines]
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