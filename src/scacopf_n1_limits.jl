function set_n1_limits!(opfdata::OPFData, options::Dict, adjustments::Dict, contingencies::Dict,
                        max_iter::Int=10, viol_scale::Float64=1.005, nonviol_scale::Float64=1.0)

    ## shortcuts
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen = length(generators)

    ## compute initial ratings at 0-contingency
    opfdata.lines.rateA .= 0.0
    point_0, M = get_scacopf_point(opfdata, options, adjustments, contingencies)
    @assert(M.status == :Optimal)
    flowmag2s_0 = get_flowmag2s(point_0, opfdata, options, false)
    ratings_0   = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
    ratings     = deepcopy(ratings_0)
    opfdata.lines.rateA .= max.(ratings, 1.0)

    ## compute initial ratings at c_id-contingency
    for c_id in keys(contingencies)
        opfd = deepcopy(opfdata)
        removed_line = remove_line!(opfd, c_id)
        point_c     = get_contingency_point(M, c_id)
        flowmag2s_c = get_flowmag2s(point_c, opfd, options, true)
        ratings_c   = get_ratings(flowmag2s_c.flowmag2, opfd.baseMVA)
        update_limits!(opfdata, ratings_c, viol_scale, nonviol_scale)
    end

    ## try solving scacopf problem and also update ratings
    solved = false; iter = 0
    while solved != true && iter <= max_iter

        ## compute ratings at 0-contingency
        point_0, M  = get_scacopf_point(opfdata, options, adjustments, contingencies)
        if M.status == :Optimal; break; end
        flowmag2s_0 = get_flowmag2s(point_0, opfdata, options)
        ratings_0   = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
        ratings     = deepcopy(ratings_0)
        update_limits!(opfdata, ratings, viol_scale, nonviol_scale)

        ## compute ratings at c_id-contingency
        for c_id in keys(contingencies)
            opfd = deepcopy(opfdata)
            removed_line = remove_line!(opfd, c_id)
            point_c      = get_contingency_point(M, c_id)
            flowmag2s_c  = get_flowmag2s(point_c, opfd, options).flowmag2; splice!(flowmag2s_c, c_id:c_id-1, 0.0)
            ratings_c    = get_ratings(flowmag2s_c, opfd.baseMVA)
            update_limits!(opfdata, ratings_c, viol_scale, nonviol_scale)
            reinstate_line!(opfd, c_id, removed_line)
        end

        ## update counter
        iter += 1
    end
end

function get_acopf_point(opfdata::OPFData, options::Dict, adjustments::Dict, warm_point=false)
    M = acopf_model(opfdata, options, adjustments)
    M = acopf_solve(M, opfdata, warm_point)
    point = Dict()
    point[:Pg] = deepcopy(getvalue(M.m[:Pg]))
    point[:Qg] = deepcopy(getvalue(M.m[:Qg]))
    point[:Vm] = deepcopy(getvalue(M.m[:Vm]))
    point[:Va] = deepcopy(getvalue(M.m[:Va]))
    return point, M
end

function get_scacopf_point(opfdata::OPFData, options::Dict, adjustments::Dict, contingencies::Dict)
    M = scacopf_model(opfdata, options, adjustments, contingencies)
    M = scacopf_solve(M, opfdata, options, contingencies)
    x_calc = MathProgBase.getsolution(M.m.internalModel)
    calc_point = Dict()
    calc_point[:Pg] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Pg)]])
    calc_point[:Qg] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Qg)]])
    calc_point[:Vm] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Vm)]])
    calc_point[:Va] = deepcopy(x_calc[[x.col for x in getindex(M.m, :Va)]])
    return calc_point, M
end

function get_contingency_point(M::OPFModel, c_id::Int)
    x_calc = MathProgBase.getsolution(M.m.internalModel)
    calc_point = Dict()
    calc_point[:Pg] = deepcopy(x_calc[[x.col for x in getindex(M.m, Symbol("Pg_$(c_id)"))]])
    calc_point[:Qg] = deepcopy(x_calc[[x.col for x in getindex(M.m, Symbol("Qg_$(c_id)"))]])
    calc_point[:Vm] = deepcopy(x_calc[[x.col for x in getindex(M.m, Symbol("Vm_$(c_id)"))]])
    calc_point[:Va] = deepcopy(x_calc[[x.col for x in getindex(M.m, Symbol("Va_$(c_id)"))]])
    return calc_point
end

function update_limits!(opfdata::OPFData, ratings::Array{Float64,1}, viol_scale::Float64, nonviol_scale::Float64=1.0, verb=false)
    viol_idx = ratings .>= opfdata.lines.rateA
    opfdata.lines.rateA[viol_idx]   .= ratings[viol_idx]   .* viol_scale
    opfdata.lines.rateA[.!viol_idx] .= max.(ratings[.!viol_idx] .* nonviol_scale, opfdata.lines.rateA[.!viol_idx])
    if verb
        println("highest limit on line $(argmax(opfdata.lines.rateA))")
    end
end
