function set_n1_limits!(opfdata::OPFData, options::Dict, feas_tol=1e-6, solve_scale=1.25, max_iter=10)
    ## shortcuts
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen = length(generators)

    ## set all line limits to zero
    opfdata.lines.rateA .= zeros(nline)

    ## set ratings for no line outages
    point_0, M_0 = get_opf_point(opfdata, options); @assert(M_0.status == :Optimal)
    flowmag2s_0  = get_flowmag2s(point_0, opfdata, options)
    ratings_0    = get_ratings(flowmag2s_0, opfdata.baseMVA)
    opfdata.lines.rateA .= ratings_0

    ## get non-islanding lines
    nonislanding_lines = get_nonislanding_lines(opfdata, options)

    ## adjust ratings for individual line outages
    point   = deepcopy(point_0)
    ratings = deepcopy(ratings_0)
    for l in nonislanding_lines
        rl = remove_line!(opfdata, l)
        adjust_feas_ratings!(opfdata, options, point, feas_tol, max_iter)
        solved, M, point = adjust_solv_ratings!(opfdata, options, point, solve_scale)
        reinstate_line!(opfdata, l, rl)
    end

    ## display
    for l in eachindex(lines)
        println("Line   From Bus    To Bus    Rating (Orig)    Rating (New)")
        @printf("%3d      %3d      %3d        %5.3f          %5.3f\n", l, lines.from[l], lines.to[l], ratings_0[l], lines.rateA[l])
    end
    return opfdata
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

function adjust_feas_ratings!(opfdata::OPFData, options::Dict, point::Dict, tol=1e-6)
    """
    modify `opfdata.lines.rateA` so that `point` is feasibile
    """
    feas, infeas_line_idx = check_feasibility(point, opfdata, options, tol)
    if feas == false
        flowmag2s   = get_flowmag2s(point, opfdata, options)
        ratings     = get_ratings(flowmag2s, opfdata.baseMVA)
        adj_ratings = max.(ratings, opfdata.lines.rateA)
        opfdata.lines.rateA .= adj_ratings
    end
    feas, infeas_line_idx = check_feasibility(point, opfdata, options, tol)
    @assert(feas == true)
    nothing
end
function adjust_solv_ratings!(opfdata::OPFData, options::Dict, point::Dict, scale=1.5, max_iter=10)
    """
    modify `opfdata.lines.rateA` so that `acopf_solve` reaches an optimal point beginning from `point`
    """
    solved, M, _ = check_solvability(point, opfdata, options)
    iter = 1
    while solved == false && iter <= max_iter
        flowmag2s   = get_flowmag2s(point, opfdata, options)
        ratings     = get_ratings(flowmag2s, opfdata.baseMVA)
        adj_ratings = max.(ratings, opfdata.lines.rateA)
        opfdata.lines.rateA .= scale .* adj_ratings
        solved, M, _ = check_solvability(point, opfdata, options)
        iter += 1
        if iter == max_iter
            @warn "Maximum number of uniform line rating scalings reached."
        end
    end
    return solved, M, point
end

function check_feasibility(M::OPFModel, opfdata::OPFData, options::Dict, tol=1e-6)
    M_eval = setup(M.m)
    x = [getvalue(M.m[:Pg]); getvalue(M.m[:Qg]); getvalue(M.m[:Vm]); getvalue(M.m[:Va])]
    c = zeros(MathProgBase.numconstr(M.m))
    c!(c, x, model=M_eval)
    nflowlim = length(opfdata.lines.rateA .> 0)
    feas = true
    c_vltgs = c[1:(end-nflowlim)]
    feas *= prod(c_vltgs .< 0)
    @assert(feas == true)
    c_flows = c[(end-nflowlim+1):end]
    feas *= prod(c .< tol)
    if feas == true
        return true, []
    else
        return false, findall(c .< tol)
    end
end

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
    VM = point[:Vm]
    VA = point[:Va]
    Y = computeAdmittanceMatrix(opfdata, options)
    Vm_lo = opfdata.buses.Vmin
    Vm_hi = opfdata.buses.Vmax
    Va_lo = -pi
    Va_hi = pi
    flowmag2s_max = get_flowmag2s(VM, VA, Y, opfdata, options)
    feas = true
    feas *= prod(VM .<= Vm_hi); feas *= prod(VM .>= Vm_lo)
    feas *= prod(VA .<= Va_hi); feas *= prod(VA .>= Va_lo)
    @assert(feas == true)
    if feas == true
        return true, []
    else
        return false, findall(c .< tol)
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
#
# function get_new_ratings(point::Dict, opfdata::OPFData, options::Dict)
#     feas, lineidx = check_feasibility(point, opfdata, options)
#     if feas == false
#         point, M  = get_opf_point(opfdata, options)
#         flowmag2s = get_flowmag2s(M, opfdata, options)
#         ratings   = get_ratings(flowmag2s, opfdata.baseMVA)
#         return ratings, flowmag2s, point, M
#     end
# end

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
