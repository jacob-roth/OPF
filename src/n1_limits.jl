function get_n1_limits(case::String, path::String,
                       initial_loading::Float64=0.46, loading_inc::Float64=0.01,
                       ratings_buffer::Float64=1.25)
    ## load data
    opfdata0 = load_case(case, path, other=false);
    nbus     = length(opfdata0.buses);
    ngen     = length(opfdata0.generators);
    nload    = sum(opfdata0.buses.bustype .== 1);
    options  = DefaultOptions();
    options[:current_rating] = true
    adjustments = DefaultAdjustments()
    max_iter = Int(ceil(initial_loading/loading_inc))

    #
    # find max loading and limits
    #
    feas = false
    loading = initial_loading
    iter = 1
    new_ratings = zeros(length(opfdata0.lines))
    while !feas
        println("==============================================")
        println("==========        ITERATION $iter       ==========")
        println("==============================================")
        ## copy original data
        opfdata = deepcopy(opfdata0)

        ## compute maximum ratings
        println("--> maximizing rateA")
        update_ratings_max!(opfdata, options)
        max_ratings = deepcopy(opfdata.lines.rateA)

        ## update loading
        println("--> modifying P & Q loading to $(round(loading*100,digits=2))%")
        update_loadings!(opfdata, options, loading, DefaultAdjPF())

        ## get contingencies
        println("--> building contingencies")
        contingencies = get_contingencies(opfdata, options)
        c = Dict()
        for k in collect(keys(contingencies))
            c[k] = contingencies[k]
        end

        c_opfdata = OPFData[]
        for k in collect(keys(contingencies))
            cc = contingencies[k]
            rl = remove_line!(opfdata, cc.asset.id);
            push!(c_opfdata, deepcopy(opfdata))
            reinstate_line!(opfdata, rl.id, rl);
        end

        ## solve scacopf
        println("--> solving max-rate sc-acopf")
        scm = scacopf_model(opfdata, options, adjustments, c)
        @time scm = scacopf_solve(scm, opfdata, options, c)
        if scm.status == :Infeasible
            loading -= loading_inc
            println("--> decreasing P & Q loading parameter to $(round(loading*100,digits=2))%")
            iter += 1
            continue
        end

        ## get dispatch point
        x = MathProgBase.getsolution(scm.m.internalModel)
        dp = Dict()
        dp[:Pg] = deepcopy(x[[x.col for x in getindex(scm.m, :Pg)]])
        dp[:Qg] = deepcopy(x[[x.col for x in getindex(scm.m, :Qg)]])
        dp[:Vm] = deepcopy(x[[x.col for x in getindex(scm.m, :Vm)]])
        dp[:Va] = deepcopy(x[[x.col for x in getindex(scm.m, :Va)]])

        ## get operating points & flow limits
        println("--> getting operating points")
        ops          = get_operating_points(scm, c)
        push!(ops, dp)
        tf           = [[true for i in eachindex(collect(keys(c)))]; false]
        flows        = [get_flowmag2s(op, c_data, options, ctf) for (op, c_data, ctf) in zip(ops, [c_opfdata; opfdata0], tf)]
        base_ratings = [get_ratings(x.flowmag2, opfdata.baseMVA) for x in flows]
        base_ratings = cat(base_ratings..., dims=2)
        adj_ratings  = map(maximum, eachrow(base_ratings))

        ## adjust ratings
        println("--> adjusting ratings")
        new_ratings = adj_ratings .* ratings_buffer
        opfdata.lines.rateA .= new_ratings
        println("--> solving adj-ratings sc-acopf")
        scm_new = scacopf_model(opfdata, options, adjustments, c)
        @time scm_new = scacopf_solve(scm, opfdata, options, c)

        ## update infeasibility
        feas = scm_new.status == :Optimal
        iter += 1
        if iter >= max_iter
            println("--> maximum iterations reached")
            break
        end
    end
    return opfdata
end