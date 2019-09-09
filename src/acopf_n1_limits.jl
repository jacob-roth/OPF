function set_n1_limits!(opfdata::OPFData, options::Dict, feas_tol=1e-6, solve_scale=1.25, max_iter=10)
    ## shortcuts
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen = length(generators)

    ## set all line limits to zero
    opfdata.lines.rateA .= zeros(nline)

    ## get initial IPOPT point (Va = slack, Vm ≈ 1, Pg & Qg from MPC casefile)
    point_init = acopf_initialPt_IPOPT_point(opfdata)

    ## set ratings for no line outages
    point_0, M_0 = get_opf_point(opfdata, options); @assert(M_0.status == :Optimal)
    flowmag2s_0  = get_flowmag2s(point_0, opfdata, options)
    ratings_0    = get_ratings(flowmag2s_0, opfdata.baseMVA)
    opfdata.lines.rateA .= ratings_0

    ## ensure N-0 solvability
    solved, M, point = adjust_solv_ratings!(opfdata, options, point_init, solve_scale, max_iter)

    ## get non-islanding lines
    nonislanding_lines = get_nonislanding_lines(opfdata, options)

    ## adjust ratings for individual line outages
    point   = deepcopy(point_0)
    ratings = deepcopy(ratings_0)
    for l in nonislanding_lines
        println("========== REMOVING LINE $l        ==========")
        rl = remove_line!(opfdata, l)
        solved, M, point = adjust_solv_ratings!(opfdata, options, point, solve_scale, max_iter)
        reinstate_line!(opfdata, l, rl)
    end

    ## display
    for l in eachindex(lines)
        println("Line   From Bus    To Bus    Rating (Orig)    Rating (New)")
        @printf("%3d      %3d        %3d          %5.3f            %5.3f\n", l, lines.from[l], lines.to[l], ratings_0[l], opfdata.lines.rateA[l])
    end
    # return opfdata
end

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

function adjust_feas_ratings!(opfdata::OPFData, options::Dict, point::Dict, tol=1e-6, buffer=0.10)
    """
    modify `opfdata.lines.rateA` so that `point` is feasibile
    """
    feas, infeas_dict = check_feasibility(point, opfdata, options, tol)
    if feas == false
        ## flow infeasibility
        flowmag2s   = get_flowmag2s(point, opfdata, options)
        ratings     = get_ratings(flowmag2s, opfdata.baseMVA)
        adj_ratings = max.(ratings, opfdata.lines.rateA)
        opfdata.lines.rateA .= adj_ratings .* (1.0 + buffer)

        ## voltage infeasibility
        buses = unique([infeas_dict[:VM_lo];
                        infeas_dict[:VM_hi];
                        infeas_dict[:VA_lo];
                        infeas_dict[:VA_hi];])
        lines = unique([findall(opfdata.lines.to .∈ Ref(buses));
                        findall(opfdata.lines.from .∈ Ref(buses))])
        opfdata.lines.rateA[lines] .*= (1.0 + buffer)
    end
    feas, infeas_dict = check_feasibility(point, opfdata, options, tol)
    @assert(feas == true)
    nothing
end
function adjust_solv_ratings!(opfdata::OPFData, options::Dict, point::Dict, buffer=0.1, max_iter=10, pct=1.0, feas_tol=1e-6)
    """
    modify `opfdata.lines.rateA` so that `acopf_solve` reaches an optimal point beginning from `point`
    by adjusting ratings according to
        1. ensuring flow feasibility
        2. increasing ratings if bus Vm falls within `pct`% of its upper/lower limits
    """
    solved, M, _ = check_solvability(point, opfdata, options)
    iter = 1
    while solved == false && iter <= max_iter
        println("~~~~~~~~~~ inner solve iteration $iter ~~~~~~~~~~")
        ## get infeasibile point from within solver
        infeas_x = MathProgBase.getsolution(M.m.internalModel)
        infeas_point = Dict()
        infeas_point[:Pg] = infeas_x[[x.col for x in getindex(M.m, :Pg)]]
        infeas_point[:Qg] = infeas_x[[x.col for x in getindex(M.m, :Qg)]]
        infeas_point[:Vm] = infeas_x[[x.col for x in getindex(M.m, :Vm)]]
        infeas_point[:Va] = infeas_x[[x.col for x in getindex(M.m, :Va)]]

        ## adjust ratings so that infeasible point is now feasible
        adjust_feas_ratings!(opfdata, options, infeas_point, feas_tol)

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
        solved, M, _ = check_solvability(infeas_point, opfdata, options)
        solved, M, _ = check_solvability(point, opfdata, options)
        iter += 1
        if iter == max_iter
            throw(ErrorException("Maximum number of line rating scalings reached."))
            # @warn "Maximum number of uniform line rating scalings reached."
        end
        if any(opfdata.lines.rateA .> 1e7)
            throw(ErrorException("Maximum line limit increase reached."))
        end
    end
    return solved, M, point
end

# function check_feasibility(M::OPFModel, opfdata::OPFData, options::Dict, tol=1e-6)
#     M_eval = setup(M.m)
#     x = [getvalue(M.m[:Pg]); getvalue(M.m[:Qg]); getvalue(M.m[:Vm]); getvalue(M.m[:Va])]
#     c = zeros(MathProgBase.numconstr(M.m))
#     c!(c, x, model=M_eval)
#     nflowlim = length(opfdata.lines.rateA .> 0)
#     feas = true
#     c_vltgs = c[1:(end-nflowlim)]
#     feas *= prod(c_vltgs .< 0)
#     @assert(feas == true)
#     c_flows = c[(end-nflowlim+1):end]
#     feas *= prod(c .< tol)
#     if feas == true
#         return true, []
#     else
#         return false, findall(c .< tol)
#     end
# end

function check_solvability(start_point::Dict, opfdata::OPFData, options::Dict)
    M = acopf_model(opfdata, options)
    M = acopf_solve(M, opfdata, start_point)
    opt = M.status == :Optimal
    point = Dict()
    point[:Pg] = deepcopy(getvalue(M.m[:Pg]))
    point[:Qg] = deepcopy(getvalue(M.m[:Qg]))
    point[:Vm] = deepcopy(getvalue(M.m[:Vm]))
    point[:Va] = deepcopy(getvalue(M.m[:Va]))
    return opt, M, point
end

function check_feasibility(point::Dict, opfdata::OPFData, options::Dict, tol=1e-6)
    ## process
    VM = point[:Vm]
    VA = point[:Va]
    Y  = computeAdmittanceMatrix(opfdata, options)
    Vm_lo = opfdata.buses.Vmin
    Vm_hi = opfdata.buses.Vmax
    Va_lo = -pi
    Va_hi = pi
    flowmag2s = get_flowmag2s(VM, VA, Y, opfdata, options)
    flowmax   = (opfdata.lines.rateA ./ opfdata.baseMVA).^2

    ## check componentwise feasibility
    feas  = true
    VM_hi = VM .<= Vm_hi
    VM_lo = VM .>= Vm_lo
    VA_hi = VA .<= Va_hi
    VA_lo = VA .>= Va_lo
    flows = flowmag2s .<= flowmax
    feas *= prod(VM_hi); feas *= prod(VM_lo)
    feas *= prod(VA_hi); feas *= prod(VA_lo)
    feas *= prod(flows)

    ## infeasibility results
    infeas_dict = Dict()
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

function remove_line!(opfdata::OPFData, l::Int64)
    lines = [x for x in opfdata.lines]
    rl = splice!(lines, l)
    opfdata.lines = StructArray(lines)

    ## adjust FromLines & ToLines
    FromLines, ToLines = mapLinesToBuses(opfdata.buses, opfdata.lines, opfdata.BusIdx)
    opfdata.FromLines .= FromLines
    opfdata.ToLines   .= ToLines
    return rl
end
function reinstate_line!(opfdata::OPFData, l::Int64, rl::MPCCases.Line)
    lines = [x for x in opfdata.lines]
    splice!(lines, (l):(l-1), [rl])
    opfdata.lines = StructArray(lines)

    ## adjust FromLines & ToLines
    FromLines, ToLines = mapLinesToBuses(opfdata.buses, opfdata.lines, opfdata.BusIdx)
    opfdata.FromLines .= FromLines
    opfdata.ToLines   .= ToLines
    nothing
end

function get_opf_point(opfdata::OPFData, options::Dict)
    M = acopf_model(opfdata, options)
    M = acopf_solve(M, opfdata)
    point = Dict()
    point[:Pg] = deepcopy(getvalue(M.m[:Pg]))
    point[:Qg] = deepcopy(getvalue(M.m[:Qg]))
    point[:Vm] = deepcopy(getvalue(M.m[:Vm]))
    point[:Va] = deepcopy(getvalue(M.m[:Va]))
    return point, M
end

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
