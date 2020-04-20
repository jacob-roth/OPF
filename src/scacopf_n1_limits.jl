function set_n1_limits!(opfdata::OPFData, options::Dict, adjustments::Dict, max_iter::Int=10, pct::Union{Int, Float64}=0.1, verbose::Bool=false)
    ## shortcuts
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen = length(generators)

    ## set ratings for no line outages
    opfdata.lines.rateA .= 0
    point_0, M_0 = get_acopf_point(opfdata, options, adjustments)
    @assert(M_0.status == :Optimal)
    flowmag2s_0 = get_flowmag2s(point_0, opfdata, options)
    ratings_0 = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
    opfdata.lines.rateA .= ratings_0

    solved = false; point = point_0; iter = 0
    contingencies = get_all_contingencies(opfdata, options)
    while solved == false && iter <= max_iter
        flowmag2s_0 = get_flowmag2s(point, opfdata, options)
        ratings_0 = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
        opfdata.lines.rateA .= update_rating_limits(opfdata, ratings_0, pct)
        if verbose == true
            println(opfdata.lines.rateA)
        end

        for c_id in keys(contingencies)
            removed_line = remove_line!(opfdata, c_id)
            flowmag2s_l = get_flowmag2s(point, opfdata, options).flowmag2
            ratings_l = get_ratings(flowmag2s_l, opfdata.baseMVA)
            opfdata.lines.rateA .= update_rating_limits(opfdata, ratings_l, pct)
            reinstate_line!(opfdata, c_id, removed_line)
        end

        point, M = get_scacopf_point(opfdata, options, adjustments, contingencies)
        if M.status == :Optimal
            solved = true
        end

        iter += 1
    end

    ## display
    # for l in eachindex(lines)
    #     println("Line   From Bus    To Bus    Rating (Orig)    Rating (New)")
    #     @printf("%3d      %3d        %3d          %5.3f            %5.3f\n", l, lines.from[l], lines.to[l], ratings_0[l], opfdata.lines.rateA[l])
    # end
    # return opfdata
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

function get_scacopf_point(opfdata::OPFData, options::Dict, adjustments::Dict, contingencies:: Dict=Dict())
    M = scacopf_model(opfdata, options, adjustments, contingencies)
    M = scacopf_solve(M, opfdata, options, contingencies)
    point = Dict()
    point[:Pg] = deepcopy(getvalue(M.m[:Pg]))
    point[:Qg] = deepcopy(getvalue(M.m[:Qg]))
    point[:Vm] = deepcopy(getvalue(M.m[:Vm]))
    point[:Va] = deepcopy(getvalue(M.m[:Va]))
    return point, M
end

function update_rating_limits(opfdata::OPFData, ratings::Array{Float64, 1}, pct::Union{Int, Float64})
    rating_limits = opfdata.lines.rateA
    updated_rating_limits = rating_limits
    violating_idx = ratings .> rating_limits
    updated_rating_limits[violating_idx] = ratings[violating_idx] * (1 + pct)
    return updated_rating_limits
end
